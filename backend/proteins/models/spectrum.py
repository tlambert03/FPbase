from __future__ import annotations

import ast
import json
import logging
import operator
from functools import cached_property, reduce
from typing import TYPE_CHECKING, Any

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
from proteins.util.spectra import interp_linear
from references.models import Reference

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import NotRequired, TypedDict

    from proteins.models import FluorState

    class D3Dict(TypedDict):
        slug: str
        key: str
        id: int
        values: list[dict[str, float]]
        peak: int | None
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
        owner_name_paths = (
            "owner_fluor__owner_name",
            "owner_filter__name",
            "owner_light__name",
            "owner_camera__name",
        )
        qs_list = []
        for owner_path in owner_name_paths:
            for s in (
                Spectrum.objects.annotate(similarity=TrigramSimilarity(owner_path, query))
                .filter(similarity__gt=threshold)
                .order_by("-similarity", owner_path)
                .distinct("similarity", owner_path)
            ):
                qs_list.append((s.similarity, s.owner))
        if qs_list:
            max_sim = max([s[0] for s in qs_list])
            qs_list = [i[1] for i in qs_list if max_sim - i[0] < 0.05]
        return qs_list


class SpectrumData(ArrayField):
    """Legacy field class - kept for migration compatibility only."""

    def __init__(self, base_field=None, size=None, **kwargs):
        if not base_field:
            base_field = ArrayField(models.FloatField(max_length=10), size=2)
        super().__init__(base_field, size, **kwargs)


def _exactly_one_of(fields: Sequence[str]) -> models.Q:
    """Return Q condition ensuring exactly one of the given fields is set (not null)."""
    conditions = (models.Q(**{f"{f}__isnull": f != field for f in fields}) for field in fields)
    return reduce(operator.or_, conditions)


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

    # Optimized storage: Y values only (wavelengths derived from min/max_wave)
    # Y values stored as float32 binary for ~75% space savings vs legacy format
    min_wave = models.SmallIntegerField(help_text="Minimum wavelength in nm (inclusive)")
    max_wave = models.SmallIntegerField(help_text="Maximum wavelength in nm (inclusive)")
    y_values = models.BinaryField(
        help_text="Y values as float32 binary array, one value per nm from min to max wave",
    )
    scale_factor = models.FloatField(
        null=True,
        blank=True,
        help_text=(
            "Physical scaling constant. Meaning depends on subtype: "
            "EX/ABS=extinction coeff (M^-1 cm^-1), EM=quantum yield, "
            "2P=peak cross-section (GM), filters/cameras=1.0"
        ),
    )
    peak_wave = models.SmallIntegerField(
        null=True,
        blank=True,
        help_text="Wavelength of peak intensity in nm",
    )

    category = models.CharField(
        max_length=1, choices=CATEGORIES, verbose_name="Spectrum Type", db_index=True
    )
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
    owner_fluor: models.ForeignKey[FluorState | None] = models.ForeignKey(
        "FluorState",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectra",
    )
    owner_filter_id: int | None
    owner_filter: models.OneToOneField[Filter | None] = models.OneToOneField(
        "Filter",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_light_id: int | None
    owner_light: models.OneToOneField[Light | None] = models.OneToOneField(
        "Light",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_camera_id: int | None
    owner_camera: models.OneToOneField[Camera | None] = models.OneToOneField(
        "Camera",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    reference_id: int | None
    reference: models.ForeignKey[Reference | None] = models.ForeignKey(
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
            # Composite index for the most common query pattern:
            # filtering by fluorophore and status
            models.Index(fields=["owner_fluor_id", "status"], name="spectrum_fluor_status_idx"),
            # Index on status for queries that only filter by approval status
            models.Index(fields=["status"], name="spectrum_status_idx"),
            # Covering index for metadata-only queries
            models.Index(fields=["status", "category", "subtype"], name="spectrum_metadata_idx"),
            # Partial index for approved fluorophore spectra (most common query)
            models.Index(
                fields=["owner_fluor_id"],
                name="spectrum_approved_fluor_idx",
                condition=models.Q(status="approved", owner_fluor_id__isnull=False),
            ),
        ]
        constraints = [
            # Ensure exactly one owner is set
            models.CheckConstraint(
                name="spectrum_single_owner",
                violation_error_message=(
                    "Spectrum must have exactly one owner (fluor, filter, light, or camera)"
                ),
                condition=_exactly_one_of(
                    ["owner_fluor", "owner_filter", "owner_light", "owner_camera"]
                ),
            ),
            # Ensure subtype is valid for category
            models.CheckConstraint(
                name="spectrum_valid_category_subtype",
                violation_error_message="Subtype is not valid for the spectrum category",
                condition=(
                    models.Q(category="d", subtype__in=["ex", "ab", "em", "2p"])
                    | models.Q(category="p", subtype__in=["ex", "ab", "em", "2p"])
                    | models.Q(category="f", subtype__in=["sp", "bx", "bs", "bp", "bm", "lp"])
                    | models.Q(category="c", subtype="qe")
                    | models.Q(category="l", subtype="pd")
                ),
            ),
        ]

    def __str__(self) -> str:
        if self.owner_fluor:
            return f"{self.owner_fluor} {self.subtype}"
        else:
            return self.name

    def save(self, *args, **kwargs) -> None:
        # Auto-set subtype if not provided and only one valid option exists
        valid_subtypes = self.category_subtypes.get(self.category, [])
        if not self.subtype and len(valid_subtypes) == 1:
            self.subtype = valid_subtypes[0]

        # Compute peak_wave from spectrum data if not already set
        if self.peak_wave is None:
            self.peak_wave = self._compute_peak_wave()

        # Cache invalidation is handled by post_save signal in fpbase.cache_utils
        super().save(*args, **kwargs)

    def clean(self) -> None:
        # Validate subtype is appropriate for category
        valid_subtypes = self.category_subtypes.get(self.category, [])
        if self.subtype and self.subtype not in valid_subtypes:
            valid_str = ", ".join(valid_subtypes)
            cat = self.get_category_display()
            raise ValidationError(
                {"subtype": f"{cat} spectrum subtype must be one of: {valid_str}"}
            )

        # Validate y_values length matches wavelength range (1nm spacing)
        # y_values is float32 binary (4 bytes per value)
        if self.y_values and self.min_wave is not None and self.max_wave is not None:
            expected_len = self.max_wave - self.min_wave + 1
            actual_len = len(self.y_values) // 4
            if actual_len != expected_len:
                raise ValidationError(
                    f"y_values length ({actual_len}) must match wavelength range ({expected_len})"
                )

    @property
    def owner(self) -> FluorState | Filter | Light | Camera:
        owners = (self.owner_fluor, self.owner_filter, self.owner_light, self.owner_camera)
        return next(x for x in owners if x is not None)

    @property
    def name(self) -> str:
        # owner fluors can have multiple spectra, so include subtype
        if self.owner_fluor:
            return f"{self.owner_fluor.owner_name} {self.subtype}"
        return str(self.owner)

    def _compute_peak_wave(self) -> int | None:
        """Compute peak wavelength from spectrum data."""
        try:
            if self.min_wave < 300:
                return self.x[
                    self.y.index(max([i for n, i in enumerate(self.y) if self.x[n] > 300]))
                ]
            else:
                # First look for Y ≈ 1.0 (the normalized peak) to avoid false 2P peaks
                y_arr = np.asarray(self.y)
                close_to_one = np.where(np.isclose(y_arr, 1.0, atol=1e-6))[0]
                if len(close_to_one) > 0:
                    return self.x[close_to_one[0]]
                # Fall back to absolute max
                return self.x[int(np.argmax(y_arr))]
        except (ValueError, IndexError):
            return None

    @property
    def step(self) -> int:
        return 1

    def scaled_data(self, scale: float) -> list[tuple[float, float]]:
        """Return spectrum data with Y values multiplied by scale."""
        return [(w, y * scale) for w, y in zip(self.x, self.y)]

    def color(self) -> str:
        return wave_to_hex(self.peak_wave)

    def waverange(self, waverange: tuple[float, float]) -> list[tuple[int, float]]:
        """Return spectrum data within the specified wavelength range."""
        start = max(self.min_wave, int(waverange[0]))
        end = min(self.max_wave, int(waverange[1]))
        start_idx = start - self.min_wave
        end_idx = end - self.min_wave + 1
        return list(zip(range(start, end + 1), self.y[start_idx:end_idx]))

    def avg(self, waverange: tuple[float, float]) -> float:
        """Return average Y value within the specified wavelength range."""
        start = max(self.min_wave, int(waverange[0]))
        end = min(self.max_wave, int(waverange[1]))
        start_idx = start - self.min_wave
        end_idx = end - self.min_wave + 1
        return np.mean(self.y[start_idx:end_idx]).item()

    def width(self, height: float = 0.5) -> tuple[int, int] | None:
        """Return (min_wave, max_wave) where Y exceeds height threshold."""
        above = np.where(np.array(self.y) > height)[0]
        if len(above) == 0:
            return None
        return (self.min_wave + above[0], self.min_wave + above[-1])

    def d3dict(self) -> D3Dict:
        if self.category == self.CAMERA:
            color = "url(#crosshatch)"
        elif self.category == self.LIGHT:
            color = "url(#wavecolor_gradient)"
        else:
            color = self.color()
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
            "color": color,
            "area": self.subtype not in (self.LP, self.BS),
            "url": self.owner.get_absolute_url(),
            "classed": f"category-{self.category} subtype-{self.subtype}",
        }

        if self.category in (self.PROTEIN, self.DYE):
            # Map subtype to (fallback_scalar_attr, extra_key, extra_attr)
            subtype_map = {
                self.EX: ("ext_coeff", "ex_max"),
                self.ABS: ("ext_coeff", "ex_max"),
                self.EM: ("qy", "em_max"),
                self.TWOP: ("twop_peak_gm", "twop_qy"),
            }
            if self.subtype in subtype_map:
                fallback, attr = subtype_map[self.subtype]
                scalar = self.scale_factor or getattr(self.owner, fallback)
                D.update({"scalar": scalar, attr: getattr(self.owner, attr)})
        return D

    def _iter_coords(self) -> Iterator[tuple[int, float]]:
        """Generator yielding (wavelength, value) pairs."""
        for i, y in enumerate(self.y):
            yield (self.min_wave + i, y)

    def d3data(self) -> list[dict[str, float]]:
        """Return spectrum data in D3.js format."""
        return [{"x": x, "y": y} for x, y in self._iter_coords()]

    def wave_value_pairs(self) -> dict[int, float]:
        """Return spectrum data as a wavelength -> value dict."""
        return dict(self._iter_coords())

    @property
    def data(self) -> list[tuple[float, float]]:
        """Return spectrum data as [[wavelength, value], ...] pairs.

        This property provides backward compatibility with code that
        expects the legacy data format.
        """
        return list(self._iter_coords())

    @data.setter
    def data(self, value: list[list[float]] | str | None) -> None:
        """Set spectrum data from [[wavelength, value], ...] pairs.

        This setter provides backward compatibility with code that
        sets data in the legacy format. Accepts either a list or a
        JSON/Python string representation.
        """
        if value is None:
            return
        # Handle string input (e.g., from form fields)
        if isinstance(value, str):
            try:
                value = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                return
        if not value or len(value) == 0:
            return
        wavelengths = [point[0] for point in value]
        y_values = [point[1] for point in value]
        self._set_spectrum_data(wavelengths, y_values)

    @cached_property
    def x(self) -> list[int]:
        """Wavelength values in nm. Derived from min_wave/max_wave."""
        return list(range(self.min_wave, self.max_wave + 1))

    @cached_property
    def y(self) -> list[float]:
        """Y values (intensity/transmission). Decoded from binary storage."""
        return self._decode_y_values().tolist()

    def _decode_y_values(self) -> np.ndarray:
        """Decode float32 binary data to list of floats."""
        if not self.y_values:
            return []
        # Handle memoryview from database
        data = bytes(self.y_values) if isinstance(self.y_values, memoryview) else self.y_values
        return np.frombuffer(data, dtype="<f4")

    @staticmethod
    def _encode_y_values(y_list: list[float]) -> bytes:
        """Encode list of floats to float32 binary data."""
        return np.asarray(y_list, dtype="<f4").tobytes()

    def _set_spectrum_data(
        self,
        wavelengths: list[int | float],
        y_values: list[float],
        scale_factor: float | None = None,
    ) -> None:
        """Set spectrum data using the new optimized storage format.

        Parameters
        ----------
        wavelengths
            Wavelength values in nm. Will be interpolated to 1nm steps if needed.
        y_values
            Y values (intensity/transmission/cross-section).
        scale_factor
            Physical scaling constant. If provided, y_values are assumed to be
            in absolute units and will be normalized by this value.
            If None, y_values are assumed to already be normalized.
        """
        # Clear cached properties
        for prop in ("x", "y"):
            if prop in self.__dict__:
                del self.__dict__[prop]

        # Interpolate to 1nm steps if needed
        if len(wavelengths) > 1:
            steps = {wavelengths[i + 1] - wavelengths[i] for i in range(len(wavelengths) - 1)}
            if steps != {1} and steps != {1.0}:
                interp_x, interp_y = interp_linear(wavelengths, y_values)
                wavelengths = list(interp_x)
                y_values = list(interp_y)

        self.min_wave = int(min(wavelengths))
        self.max_wave = int(max(wavelengths))

        # Normalize if scale_factor provided
        if scale_factor is not None and scale_factor != 0:
            y_values = [v / scale_factor for v in y_values]
            self.scale_factor = scale_factor
        elif scale_factor is not None:
            self.scale_factor = scale_factor

        self.y_values = self._encode_y_values(y_values)

    def change_y(self, value: list[float]) -> None:
        """Change Y values in place."""
        if not isinstance(value, list):
            raise TypeError("Y values must be a python list")

        expected_len = self.max_wave - self.min_wave + 1
        if len(value) != expected_len:
            raise ValueError(
                f"Error: array length {len(value)} must match expected {expected_len}"
            )

        # Clear cached property before modifying data
        if "y" in self.__dict__:
            del self.__dict__["y"]

        self.y_values = self._encode_y_values(value)

    def spectrum_img(self, fmt="svg", **kwargs: Any):
        """Generate a static image of this spectrum using matplotlib."""

        kwargs.setdefault("xlim", (self.min_wave - 10, self.max_wave + 10))
        # Use the existing spectra_fig function with this single spectrum
        return spectra_fig([self], fmt, **kwargs)

    def get_absolute_url(self):
        return reverse("proteins:spectra") + f"?s={self.id}"


# ###########################################################################
# Spectrum Owner Classes
# ###########################################################################


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
    tavg = models.FloatField(
        blank=True, null=True, validators=[MinValueValidator(0), MaxValueValidator(1)]
    )
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
