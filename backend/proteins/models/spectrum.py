from __future__ import annotations

import ast
import json
import logging
from functools import cached_property
from typing import TYPE_CHECKING, Any

import django.forms
import numpy as np
from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.search import TrigramSimilarity
from django.core.cache import cache
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.db.models import Case, CharField, F, IntegerField, QuerySet, Value, When
from django.urls import reverse
from django.utils.text import slugify
from model_utils import Choices
from model_utils.managers import QueryManager
from model_utils.models import StatusModel, TimeStampedModel

from fpbase.cache_utils import SPECTRA_CACHE_KEY
from proteins.models.mixins import AdminURLMixin, Authorable, Product
from proteins.util.helpers import spectra_fig, wave_to_hex
from proteins.util.spectra import interp_linear, norm2one, norm2P, step_size
from references.models import Reference

if TYPE_CHECKING:
    from typing import NotRequired, TypedDict

    from proteins.models import Fluorophore

    class D3Dict(TypedDict):
        slug: str
        key: str
        id: int
        values: list[dict[str, float]]
        peak: float | bool
        minwave: float
        maxwave: float
        category: str
        type: str
        color: str
        area: bool
        url: str | None
        classed: str

        scalar: NotRequired[float | None]
        ex_max: NotRequired[float | None]
        em_max: NotRequired[float | None]
        twop_qy: NotRequired[float | None]


logger = logging.getLogger(__name__)


class SpectrumOwner(Authorable, TimeStampedModel):
    name = models.CharField(max_length=100)  # required
    slug = models.SlugField(
        max_length=128, unique=True, help_text="Unique slug for the %(class)"
    )  # calculated at save

    class Meta:
        abstract = True
        ordering = ["name"]

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.slug}>"

    def save(self, *args, **kwargs):
        self.slug = self.makeslug()
        super().save(*args, **kwargs)

    def makeslug(self):
        return slugify(self.name.replace("/", "-"))

    def d3_dicts(self) -> list[D3Dict]:
        return [spect.d3dict() for spect in self.spectra.all()]


class Camera(SpectrumOwner, Product):
    manufacturer = models.CharField(max_length=128, blank=True)

    if TYPE_CHECKING:
        spectrum: Spectrum
        microscopes: models.QuerySet
        optical_configs: models.QuerySet


class Light(SpectrumOwner, Product):
    manufacturer = models.CharField(max_length=128, blank=True)

    if TYPE_CHECKING:
        spectrum: Spectrum
        microscopes: models.QuerySet
        optical_configs: models.QuerySet


def sorted_ex2em(filterset):
    def _sort(stype):
        return Spectrum.category_subtypes[Spectrum.FILTER].index(stype)

    return sorted(filterset, key=lambda x: _sort(x.subtype))


def get_cached_spectra_info() -> str:
    """Get cached spectra info JSON, populating cache if needed.

    Returns a JSON string of spectra data. The cache is invalidated by signals
    when any related models change. Use with @condition decorator for ETags.
    """
    cached = cache.get(SPECTRA_CACHE_KEY)
    if not cached:
        cached = json.dumps({"data": {"spectra": get_spectra_list()}})
        # Cache indefinitely, rely on signals for invalidation
        if not cache.add(SPECTRA_CACHE_KEY, cached, None):
            # Another process set it first; get the value again
            cached = cache.get(SPECTRA_CACHE_KEY)
    return cached


def get_spectra_list(query_set: QuerySet | None = None, **filters: str) -> list[dict]:
    """Fetch spectra with polymorphic owner info in a single optimized query."""
    if query_set is None:
        filters.setdefault("status", Spectrum.STATUS.approved)
        qs = Spectrum.objects.filter(**filters)
    else:
        qs = query_set

    owners = ["fluor", "filter", "light", "camera"]

    # Helper to create CASE statement for polymorphic owner fields
    def owner_case(field: str, **owner_fields) -> Case:
        """Create CASE/WHEN for selecting the right owner field."""
        whens = []
        for o in owners:
            field_expr = owner_fields.get(o) if o in owner_fields else F(f"owner_{o}__{field}")
            whens.append(When(**{f"owner_{o}_id__isnull": False}, then=field_expr))
        return Case(*whens, output_field=CharField() if field != "id" else IntegerField())

    # Fetch only needed fields and annotate with owner info
    qs = (
        qs.only("id", "category", "subtype", *(f"owner_{o}_id" for o in owners))
        .select_related("owner_fluor", "owner_filter", "owner_light", "owner_camera")
        .annotate(
            owner_id=owner_case("id"),
            owner_slug=owner_case("slug"),
            # Fluorophore uses cached owner_name, others use direct name field
            owner_name=owner_case(
                "name",
                fluor=F("owner_fluor__owner_name"),
            ),
            # Fluorophore doesn't have URL, others do
            owner_url=owner_case(
                "url",
                fluor=Value(""),  # Fluorophore has no URL field
            ),
        )
        .order_by("owner_name")
        .values("id", "category", "subtype", "owner_id", "owner_slug", "owner_name", "owner_url")
    )

    return [
        {
            "id": item["id"],
            "category": item["category"],
            "subtype": item["subtype"],
            "owner": {k: item[f"owner_{k}"] for k in ["id", "slug", "name", "url"]},
        }
        for item in qs
    ]


class SpectrumManager(models.Manager):
    def get_queryset(self):
        # by default, only include approved spectra
        return super().get_queryset().filter(status=Spectrum.STATUS.approved)

    def all_objects(self):
        return super().get_queryset()

    def fluor_slugs(self):
        """Get all fluorophore (State + DyeState) slugs."""
        return (
            self.get_queryset()
            .exclude(owner_fluor=None)
            .values_list("owner_fluor__slug", "owner_fluor__owner_name")
            .distinct()
        )

    def fluorlist(self):
        """Get list of all fluorophores (States + DyeStates) with spectra."""
        vallist = [
            "category",
            "subtype",
            "owner_fluor__slug",
            "owner_fluor__owner_name",
        ]
        Q = (
            self.get_queryset()
            .filter(models.Q(category=Spectrum.DYE) | models.Q(category=Spectrum.PROTEIN))
            .exclude(owner_fluor=None)
            .values(*vallist)
            .distinct("owner_fluor__slug")
        )

        out = []
        for v in Q:
            out.append(
                {
                    "category": v["category"],
                    "subtype": v["subtype"],
                    "slug": v["owner_fluor__slug"],
                    "name": v["owner_fluor__owner_name"],
                }
            )
        return sorted(out, key=lambda k: k["name"])

    def filter_owner(self, slug: str) -> QuerySet[Spectrum]:
        qs = self.none()
        A = ("owner_fluor", "owner_filter", "owner_light", "owner_camera")
        for ownerclass in A:
            qs = qs | self.get_queryset().filter(**{ownerclass + "__slug": slug})
        return qs

    def find_similar_owners(self, query, threshold=0.4):
        A = (
            "owner_fluor__owner_name",  # Search on cached parent name (Protein/Dye)
            "owner_filter__name",
            "owner_light__name",
            "owner_camera__name",
        )
        qs_list = []
        for ownerclass in A:
            for s in (
                Spectrum.objects.annotate(similarity=TrigramSimilarity(ownerclass, query))
                .filter(similarity__gt=threshold)
                .order_by("-similarity", ownerclass)
                .distinct("similarity", ownerclass)
            ):
                qs_list.append((s.similarity, s.owner))
        if qs_list:
            max_sim = max([s[0] for s in qs_list])
            qs_list = [i[1] for i in qs_list if max_sim - i[0] < 0.05]
        return qs_list


class SpectrumData(ArrayField):
    def __init__(self, base_field=None, size=None, **kwargs):
        if not base_field:
            base_field = ArrayField(models.FloatField(max_length=10), size=2)
        super().__init__(base_field, size, **kwargs)

    def formfield(self, **kwargs):
        defaults = {
            "max_length": self.size,
            "widget": django.forms.Textarea(attrs={"cols": "102", "rows": "15"}),
            "form_class": django.forms.CharField,
        }
        defaults.update(kwargs)
        return models.Field().formfield(**defaults)

    def to_python(self, value):
        if not value:
            return None
        if isinstance(value, list):
            return value
        if isinstance(value, str):
            try:
                return ast.literal_eval(value)
            except Exception as e:
                raise ValidationError("Invalid input for spectrum data") from e

    def value_to_string(self, obj):
        return json.dumps(self.value_from_object(obj))

    def clean(self, value, model_instance):
        if not value:
            return None
        value = super().clean(value, model_instance)
        step = step_size(value)
        if step > 10 and len(value) < 10:
            raise ValidationError("insufficient data")
        if step != 1:
            try:
                # TODO:  better choice of interpolation
                value = [list(i) for i in zip(*interp_linear(*zip(*value)))]
            except ValueError as e:
                raise ValidationError(f"could not properly interpolate data: {e}") from e
        return value

    def validate(self, value, model_instance):
        super().validate(value, model_instance)
        for elem in value:
            if len(elem) != 2:
                raise ValidationError("All elements in Spectrum list must have two items")
            if not all(isinstance(n, int | float) for n in elem):
                raise ValidationError("All items in Spectrum list elements must be numbers")


class Spectrum(Authorable, StatusModel, TimeStampedModel, AdminURLMixin):
    STATUS = Choices("approved", "pending", "rejected")

    DYE = "d"
    PROTEIN = "p"
    LIGHT = "l"
    FILTER = "f"
    CAMERA = "c"
    CATEGORIES = (
        (DYE, "Dye"),
        (PROTEIN, "Protein"),
        (LIGHT, "Light Source"),
        (FILTER, "Filter"),
        (CAMERA, "Camera"),
    )

    EX = "ex"
    ABS = "ab"
    EM = "em"
    TWOP = "2p"
    BP = "bp"
    BPX = "bx"  # bandpass excitation filter
    BPM = "bm"
    SP = "sp"
    LP = "lp"
    BS = "bs"  # dichroic
    QE = "qe"
    PD = "pd"
    SUBTYPE_CHOICES = (
        (EX, "Excitation"),  # for fluorophores
        (ABS, "Absorption"),  # for fluorophores
        (EM, "Emission"),  # for fluorophores
        (TWOP, "Two Photon Abs"),  # for fluorophores
        (BP, "Bandpass"),  # only for filters
        (BPX, "Bandpass-Ex"),  # only for filters
        (BPM, "Bandpass-Em"),  # only for filters
        (SP, "Shortpass"),  # only for filters
        (LP, "Longpass"),  # only for filters
        (BS, "Beamsplitter"),  # only for filters
        (QE, "Quantum Efficiency"),  # only for cameras
        (PD, "Power Distribution"),  # only for light sources
    )

    # not all subtypes are valid for all categories
    category_subtypes = {
        DYE: [EX, ABS, EM, TWOP],
        PROTEIN: [EX, ABS, EM, TWOP],
        FILTER: [SP, BPX, BS, BP, BPM, LP],
        CAMERA: [QE],
        LIGHT: [PD],
    }

    data = SpectrumData()
    category = models.CharField(max_length=1, choices=CATEGORIES, verbose_name="Spectrum Type", db_index=True)
    subtype = models.CharField(
        max_length=2,
        choices=SUBTYPE_CHOICES,
        verbose_name="Spectrum Subtype",
        db_index=True,
    )
    ph = models.FloatField(null=True, blank=True, verbose_name="pH")  # pH of measurement
    solvent = models.CharField(max_length=128, blank=True)

    # I was swayed to avoid Generic Foreign Keys by this article
    # https://lukeplant.me.uk/blog/posts/avoid-django-genericforeignkey/
    # Fluorophore encompasses both State (ProteinState) and DyeState via MTI
    owner_fluor_id: int | None
    owner_fluor = models.ForeignKey["Fluorophore | None"](
        "Fluorophore",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectra",
    )
    owner_filter_id: int | None
    owner_filter = models.OneToOneField["Filter | None"](
        "Filter",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_light_id: int | None
    owner_light = models.OneToOneField["Light | None"](
        "Light",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_camera_id: int | None
    owner_camera = models.OneToOneField["Camera | None"](
        "Camera",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    reference_id: int | None
    reference = models.ForeignKey["Reference | None"](
        Reference,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="spectra",
    )
    source = models.CharField(max_length=128, blank=True, help_text="Source of the spectra data")

    objects: SpectrumManager = SpectrumManager()
    fluorophores = QueryManager(models.Q(category=DYE) | models.Q(category=PROTEIN))
    proteins = QueryManager(category=PROTEIN)
    dyes = QueryManager(category=DYE)
    lights = QueryManager(category=LIGHT)
    filters = QueryManager(category=FILTER)
    cameras = QueryManager(category=CAMERA)

    class Meta:
        verbose_name_plural = "spectra"
        indexes = [
            # Composite index for the most common query pattern: filtering by fluorophore and status
            models.Index(fields=["owner_fluor_id", "status"], name="spectrum_fluor_status_idx"),
            # Index on status for queries that only filter by approval status
            models.Index(fields=["status"], name="spectrum_status_idx"),
        ]

    def __str__(self):
        if self.owner_fluor:
            return f"{self.owner_fluor} {self.subtype}"
        else:
            return self.name

    def save(self, *args, **kwargs):
        # FIXME: figure out why self.full_clean() throws validation error with
        # 'data cannot be null' ... even if data is provided...
        self.full_clean()
        if not any(self.owner_set):
            raise ValidationError("Spectrum must have an owner!")
        if sum(bool(x) for x in self.owner_set) > 1:
            raise ValidationError("Spectrum must have only one owner!")
        # self.category = self.owner.__class__.__name__.lower()[0]
        # Cache invalidation is handled by post_save signal in fpbase.cache_utils
        super().save(*args, **kwargs)

    def _norm2one(self):
        try:
            if self.subtype == self.TWOP:
                y, self._peakval2p, maxi = norm2P(self.y)
                self._peakwave2p = self.x[maxi]
                self.change_y(y)
            else:
                self.change_y(norm2one(self.y))
        except Exception:
            logger.exception("Error normalizing spectrum data")

    def clean(self):
        # model-wide validation after individual fields have been cleaned
        errors = {}
        if self.category == self.CAMERA:
            self.subtype = self.QE
        if self.category == self.LIGHT:
            self.subtype = self.PD
        if self.category in self.category_subtypes:
            if self.subtype not in self.category_subtypes[self.category]:
                errors.update(
                    {
                        "subtype": "{} spectrum subtype must be{} {}".format(
                            self.get_category_display(),
                            "" if len(self.category_subtypes[self.category]) > 1 else "  one of:",
                            " ".join(self.category_subtypes[self.category]),
                        )
                    }
                )

        if errors:
            raise ValidationError(errors)

        if self.data:
            if self.category == self.PROTEIN:
                self._norm2one()
            elif (max(self.y) > 1.5) or (max(self.y) < 0.1):
                if self.category in (self.FILTER, self.CAMERA) and (10 < max(self.y) < 101):
                    # assume 100% scale
                    self.change_y([round(yy / 100, 4) for yy in self.y])
                else:
                    self._norm2one()

        if errors:
            raise ValidationError(errors)

    @property
    def owner_set(self) -> list[Fluorophore | Filter | Light | Camera | None]:
        return [
            self.owner_fluor,
            self.owner_filter,
            self.owner_light,
            self.owner_camera,
        ]

    @property
    def owner(self) -> Fluorophore | Filter | Light | Camera | None:
        return next((x for x in self.owner_set if x), None)
        #  raise AssertionError("No owner is set")

    @property
    def name(self) -> str:
        # this method allows the protein name to have changed in the meantime
        if self.owner_fluor:
            # Check if it's a State (ProteinState) by checking for protein attribute
            if hasattr(self.owner_fluor, "protein"):
                if self.owner_fluor.name == "default":
                    return f"{self.owner_fluor.protein} {self.subtype}"
                else:
                    return f"{self.owner_fluor} {self.subtype}"
            else:
                # It's a DyeState
                return f"{self.owner} {self.subtype}"
        elif self.owner_filter:
            return str(self.owner)
        else:
            return str(self.owner)

    @property
    def peak_wave(self):
        try:
            if self.min_wave < 300:
                return self.x[self.y.index(max([i for n, i in enumerate(self.y) if self.x[n] > 300]))]
            else:
                try:
                    # first look for the value 1
                    # this is to avoid false 2P peaks
                    return self.x[self.y.index(1)]
                except ValueError:
                    return self.x[self.y.index(max(self.y))]
        except ValueError:
            return False

    @property
    def min_wave(self):
        return self.data[0][0]

    @property
    def max_wave(self):
        return self.data[-1][0]

    @property
    def step(self):
        s = set()
        for i in range(len(self.x) - 1):
            s.add(self.x[i + 1] - self.x[i])
        if len(s) > 1:  # multiple step sizes
            return False
        return next(iter(s))

    def scaled_data(self, scale):
        return [[n[0], n[1] * scale] for n in self.data]

    def color(self):
        return wave_to_hex(self.peak_wave)

    def waverange(self, waverange):
        assert len(waverange) == 2, "waverange argument must be an iterable of len 2"
        return [d for d in self.data if waverange[0] <= d[0] <= waverange[1]]

    def avg(self, waverange):
        d = self.waverange(waverange)
        return np.mean([i[1] for i in d])

    def width(self, height=0.5) -> tuple[float, float] | None:
        try:
            upindex = next(x[0] for x in enumerate(self.y) if x[1] > height)
            downindex = len(self.y) - next(x[0] for x in enumerate(reversed(self.y)) if x[1] > height)
            return (self.x[upindex], self.x[downindex])
        except Exception:
            return None

    def d3dict(self) -> D3Dict:
        D: D3Dict = {
            "slug": self.owner.slug,
            "key": self.name,
            "id": self.id,
            "values": self.d3data(),
            "peak": self.peak_wave,
            "minwave": self.min_wave,
            "maxwave": self.max_wave,
            "category": self.category,
            "type": self.subtype if self.subtype != self.ABS else self.EX,
            "color": self.color(),
            "area": False if self.subtype in (self.LP, self.BS) else True,
            "url": self.owner.get_absolute_url(),
            "classed": f"category-{self.category} subtype-{self.subtype}",
        }

        if self.category == self.CAMERA:
            D["color"] = "url(#crosshatch)"
        elif self.category == self.LIGHT:
            D["color"] = "url(#wavecolor_gradient)"

        if self.category in (self.PROTEIN, self.DYE):
            if self.subtype == self.EX or self.subtype == self.ABS:
                D.update({"scalar": self.owner.ext_coeff, "ex_max": self.owner.ex_max})
            elif self.subtype == self.EM:
                D.update({"scalar": self.owner.qy, "em_max": self.owner.em_max})
            elif self.subtype == self.TWOP:
                D.update({"scalar": self.owner.twop_peak_gm, "twop_qy": self.owner.twop_qy})
        return D

    def d3data(self):
        return [{"x": elem[0], "y": elem[1]} for elem in self.data]

    def wave_value_pairs(self):
        output = {}
        # arrayLength = len(self.data)
        for elem in self.data:
            output[elem[0]] = elem[1]
        return output

    @cached_property
    def x(self):
        """Extract x values from data. Cached to avoid repeated allocations."""
        return [i[0] for i in self.data]

    @cached_property
    def y(self):
        """Extract y values from data. Cached to avoid repeated allocations."""
        return [i[1] for i in self.data]

    def change_x(self, value):
        if not isinstance(value, list):
            raise TypeError("X values must be a python list")
        if len(value) != len(self.data):
            raise ValueError("Error: array length must match existing data")
        # Clear cached property before modifying data
        if "x" in self.__dict__:
            del self.__dict__["x"]
        for i in range(len(value)):
            self.data[i][0] = value[i]

    def change_y(self, value):
        if not isinstance(value, list):
            raise TypeError("Y values must be a python list")
        if len(value) != len(self.data):
            raise ValueError("Error: array length must match existing data")
        # Clear cached property before modifying data
        if "y" in self.__dict__:
            del self.__dict__["y"]
        for i in range(len(value)):
            self.data[i][1] = value[i]

    def spectrum_img(self, fmt="svg", **kwargs: Any):
        """Generate a static image of this spectrum using matplotlib."""

        kwargs.setdefault("xlim", (self.min_wave - 10, self.max_wave + 10))
        # Use the existing spectra_fig function with this single spectrum
        return spectra_fig([self], fmt, **kwargs)

    def get_absolute_url(self):
        return reverse("proteins:spectra") + f"?s={self.id}"


class Filter(SpectrumOwner, Product):
    # copied for convenience
    BP = Spectrum.BP
    BPX = Spectrum.BPX
    BPM = Spectrum.BPM
    SP = Spectrum.SP
    LP = Spectrum.LP
    BS = Spectrum.BS

    bandcenter = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(200), MaxValueValidator(1600)],
    )
    bandwidth = models.PositiveSmallIntegerField(
        blank=True, null=True, validators=[MinValueValidator(0), MaxValueValidator(900)]
    )
    edge = models.FloatField(
        null=True,
        blank=True,
        validators=[MinValueValidator(300), MaxValueValidator(1600)],
    )
    tavg = models.FloatField(blank=True, null=True, validators=[MinValueValidator(0), MaxValueValidator(1)])
    aoi = models.PositiveSmallIntegerField(
        blank=True, null=True, validators=[MinValueValidator(0), MaxValueValidator(90)]
    )

    if TYPE_CHECKING:
        from proteins.models import OpticalConfig

        spectrum: Spectrum
        optical_configs: models.QuerySet[OpticalConfig]

    class Meta:
        ordering = ["bandcenter"]

    @property
    def subtype(self):
        return self.spectrum.subtype

    @subtype.setter
    def subtype(self, x):
        self.spectrum.subtype = x
        self.spectrum.save()

    def save(self, *args, **kwargs):
        if "/" in self.name:
            try:
                if self.name.count("/") == 1:
                    w = self.name.split("/")[0]
                    self.bandcenter = int("".join([i for i in w[-4:] if i.isdigit()]))
                    w = self.name.split("/")[1].split(" ")[0]
                    self.bandwidth = int("".join([i for i in w[:4] if i.isdigit()]))
            except Exception:
                pass

        if self.bandcenter and self.bandwidth:
            try:
                wrange = (
                    (self.bandcenter - self.bandwidth / 2) + 2,
                    (self.bandcenter + self.bandwidth / 2) - 2,
                )
                self.tavg = self.spectrum.avg(wrange)
            except Exception:
                pass
        super().save(*args, **kwargs)
