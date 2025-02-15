import ast
import json
import logging

import numpy as np
from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.search import TrigramSimilarity
from django.core.cache import cache
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.forms import CharField, Textarea
from django.urls import reverse
from django.utils.text import slugify
from model_utils import Choices
from model_utils.managers import QueryManager
from model_utils.models import StatusModel, TimeStampedModel

from references.models import Reference

from ..util.helpers import wave_to_hex
from ..util.spectra import interp_linear, interp_univar, norm2one, norm2P, step_size
from .mixins import AdminURLMixin, Authorable, Product

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

    def d3_dicts(self):
        return [spect.d3dict() for spect in self.spectra.all()]


class Camera(SpectrumOwner, Product):
    manufacturer = models.CharField(max_length=128, blank=True)


class Light(SpectrumOwner, Product):
    manufacturer = models.CharField(max_length=128, blank=True)


def sorted_ex2em(filterset):
    def _sort(stype):
        return Spectrum.category_subtypes[Spectrum.FILTER].index(stype)

    return sorted(filterset, key=lambda x: _sort(x.subtype))


SPECTRA_CACHE_KEY = "spectra_sluglist"


def get_cached_spectra_info(timeout=60 * 60):
    spectrainfo = cache.get(SPECTRA_CACHE_KEY)
    if not spectrainfo:
        spectrainfo = json.dumps({"data": {"spectra": Spectrum.objects.sluglist()}})
        cache.set(SPECTRA_CACHE_KEY, spectrainfo, timeout)
    return spectrainfo


class SpectrumManager(models.Manager):
    def get_queryset(self):
        # by default, only include approved spectra
        return super().get_queryset().filter(status=Spectrum.STATUS.approved)

    def all_objects(self):
        return super().get_queryset()

    def state_slugs(self):
        L = (
            self.get_queryset()
            .exclude(owner_state=None)
            .values_list("owner_state__slug", "owner_state__protein__name", "owner_state__name")
            .distinct()
        )
        return [(slug, prot if state == "default" else f"{prot} ({state})") for slug, prot, state in L]

    def dye_slugs(self):
        return (
            self.get_queryset().filter(category=self.DYE).values_list("owner_dye__slug", "owner_dye__name").distinct()
        )

    def sluglist(self, filters: dict | None = None):
        """probably using this one going forward for spectra page"""

        owners = ["state", "dye", "filter", "light", "camera"]
        vals = [
            "id",
            "category",
            "subtype",
            "owner_state__protein__name",
            "owner_dye__url",
            "owner_filter__url",
            "owner_camera__url",
            "owner_light__url",
            "owner_state__protein__slug",
        ]
        for suffix in ["slug", "id", "name"]:
            for owner in owners:
                vals.append(f"owner_{owner}__{suffix}")
        Q = self.get_queryset()
        if filters:
            Q = Q.filter(**filters)
        Q = Q.values(*vals)

        out = []
        for v in Q:
            slug = (
                v["owner_state__slug"]
                or v["owner_dye__slug"]
                or v["owner_filter__slug"]
                or v["owner_light__slug"]
                or v["owner_camera__slug"]
            )
            name = (
                v["owner_dye__name"]
                or v["owner_filter__name"]
                or v["owner_light__name"]
                or v["owner_camera__name"]
                or None
            )
            url = (
                v["owner_dye__url"]
                or v["owner_filter__url"]
                or v["owner_state__protein__slug"]
                or v["owner_light__url"]
                or v["owner_camera__url"]
                or None
            )
            owner_id = (
                v["owner_state__id"]
                or v["owner_dye__id"]
                or v["owner_filter__id"]
                or v["owner_light__id"]
                or v["owner_camera__id"]
                or None
            )
            if not name:
                prot = v["owner_state__protein__name"]
                state = v["owner_state__name"]
                name = prot if state == "default" else f"{prot} ({state})"
            out.append(
                {
                    "id": v["id"],
                    "category": v["category"],
                    "subtype": v["subtype"],
                    "owner": {"slug": slug, "name": name, "id": owner_id, "url": url},
                }
            )
        return sorted(out, key=lambda k: k["owner"]["name"])

    # FIXME:  Stupid dumb dumb
    def fluorlist(self, withdyes=True):
        vallist = [
            "category",
            "subtype",
            "owner_state__protein__name",
            "owner_state__slug",
            "owner_state__name",
        ]
        distinct = ["owner_state__slug"]
        if withdyes:
            vallist += ["owner_dye__slug", "owner_dye__name"]
            distinct += ["owner_dye__slug"]
        Q = (
            self.get_queryset()
            .filter(models.Q(category=Spectrum.DYE) | models.Q(category=Spectrum.PROTEIN))
            .values(*vallist)
            .distinct(*distinct)
        )

        out = []
        for v in Q:
            slug = v.get("owner_state__slug") or v.get("owner_dye__slug")
            name = v.get("owner_dye__name", None)
            if not name:
                prot = v["owner_state__protein__name"]
                state = v["owner_state__name"]
                name = prot if state == "default" else f"{prot} ({state})"
            out.append(
                {
                    "category": v["category"],
                    "subtype": v["subtype"],
                    "slug": slug,
                    "name": name,
                }
            )
        return sorted(out, key=lambda k: k["name"])

    # def owner_slugs(self):
    #     ''' unused? '''
    #     L = self.get_queryset().exclude(owner_state=None).values_list(
    #         'owner_state__slug', 'owner_state__protein__name',
    #         'owner_state__name', 'category', 'subtype').distinct()
    #     out = [(slug, prot, cat if state == 'default' else '{} ({})'.format(prot, state))
    #            for slug, prot, state, cat in L]

    #     for n in ('dye', 'light', 'filter', 'camera'):
    #         out += self.get_queryset().exclude(**{'owner_' + n: None}).values_list(
    #             *['owner_' + n + '__slug', 'owner_' + n + '__name', 'category', 'subtype']).distinct()

    #     out.sort(key=lambda x: x[1], reverse=False)

    #     return out

    def filter_owner(self, slug):
        qs = self.none()
        A = ("owner_state", "owner_dye", "owner_filter", "owner_light", "owner_camera")
        for ownerclass in A:
            qs = qs | self.get_queryset().filter(**{ownerclass + "__slug": slug})
        return qs

    def find_similar_owners(self, query, threshold=0.4):
        A = (
            "owner_state__protein__name",
            "owner_dye__name",
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
            "widget": Textarea(attrs={"cols": "102", "rows": "15"}),
            "form_class": CharField,
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

    def clean(self, raw_value, model_instance):
        if not raw_value:
            return None
        raw_value = super().clean(raw_value, model_instance)
        step = step_size(raw_value)
        if step > 10 and len(raw_value) < 10:
            raise ValidationError("insufficient data")
        if step != 1:
            try:
                # TODO:  better choice of interpolation
                raw_value = [list(i) for i in zip(*interp_linear(*zip(*raw_value)))]
            except ValueError as e:
                raise ValidationError(f"could not properly interpolate data: {e}") from e
        return raw_value

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
    owner_state = models.ForeignKey("State", null=True, blank=True, on_delete=models.CASCADE, related_name="spectra")
    owner_dye = models.ForeignKey("Dye", null=True, blank=True, on_delete=models.CASCADE, related_name="spectra")
    owner_filter = models.OneToOneField(
        "Filter",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_light = models.OneToOneField(
        "Light",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    owner_camera = models.OneToOneField(
        "Camera",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="spectrum",
    )
    reference = models.ForeignKey(
        Reference,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="spectra",
    )
    source = models.CharField(max_length=128, blank=True, help_text="Source of the spectra data")

    objects = SpectrumManager()
    fluorophores = QueryManager(models.Q(category=DYE) | models.Q(category=PROTEIN))
    proteins = QueryManager(category=PROTEIN)
    dyes = QueryManager(category=DYE)
    lights = QueryManager(category=LIGHT)
    filters = QueryManager(category=FILTER)
    cameras = QueryManager(category=CAMERA)

    class Meta:
        verbose_name_plural = "spectra"

    def __str__(self):
        if self.owner_state:
            return "{} {}".format(self.owner_state if self.owner_state else "unowned", self.subtype)
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
        cache.delete(SPECTRA_CACHE_KEY)
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

    def _interpolated_data(self, method=None, **kwargs):
        if not method or method.lower() == "linear":
            return [list(i) for i in zip(*interp_linear(*zip(*self.data)))]
        elif method == "univar":
            return [list(i) for i in zip(*interp_univar(*zip(*self.data), **kwargs))]

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
    def owner_set(self):
        return [
            self.owner_state,
            self.owner_dye,
            self.owner_filter,
            self.owner_light,
            self.owner_camera,
        ]

    @property
    def owner(self):
        return next((x for x in self.owner_set if x), None)
        #  raise AssertionError("No owner is set")

    @property
    def name(self):
        # this method allows the protein name to have changed in the meantime
        if self.owner_state:
            if self.owner_state.name == "default":
                return f"{self.owner_state.protein} {self.subtype}"
            else:
                return f"{self.owner_state} {self.subtype}"
        elif self.owner_dye:
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

    def width(self, height=0.5):
        try:
            upindex = next(x[0] for x in enumerate(self.y) if x[1] > height)
            downindex = len(self.y) - next(x[0] for x in enumerate(reversed(self.y)) if x[1] > height)
            return (self.x[upindex], self.x[downindex])
        except Exception:
            return False

    def d3dict(self):
        D = {
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
                D.update({"scalar": self.owner.twop_peakGM, "twop_qy": self.owner.twop_qy})
        return D

    def d3data(self):
        return [{"x": elem[0], "y": elem[1]} for elem in self.data]
        # out = [{'x': w, 'y': 0} for w in range(300, int(self.min_wave))]
        # out.extend([{'x': elem[0], 'y': elem[1]} for elem in self.data])
        # out.extend([{'x': w, 'y': 0} for w in range(int(self.max_wave + 1), 801)])
        # return out

    def wave_value_pairs(self):
        output = {}
        # arrayLength = len(self.data)
        for elem in self.data:
            output[elem[0]] = elem[1]
        return output

    @property
    def x(self):
        self._x = []
        for i in self.data:
            self._x.append(i[0])
        return self._x

    @property
    def y(self):
        self._y = []
        for i in self.data:
            self._y.append(i[1])
        return self._y

    # def change_x(self, value):
    #     if not isinstance(value, list):
    #         raise TypeError("X values must be a python list")
    #     if len(value) != len(self.data):
    #         raise ValueError("Error: array length must match existing data")
    #     for i in range(len(value)):
    #         self.data[i][0] = value[i]

    def change_y(self, value):
        if not isinstance(value, list):
            raise TypeError("Y values must be a python list")
        if len(value) != len(self.data):
            raise ValueError("Error: array length must match existing data")
        for i in range(len(value)):
            self.data[i][1] = value[i]

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
