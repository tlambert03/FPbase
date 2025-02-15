import json
import urllib.parse

from django.contrib.postgres.fields import ArrayField
from django.core.cache import cache
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.urls import reverse
from django.utils.functional import cached_property

from ..util.efficiency import spectral_product
from ..util.helpers import shortuuid
from .collection import OwnedCollection
from .spectrum import Camera, Filter, Light, sorted_ex2em


class Microscope(OwnedCollection):
    """A microscope or other collection of optical configurations

    Items are intended to be added as configs, but can be directly added as
    filters/cameras/lights etc as well.

        configs: optical configs in the collection
                 (all spectra will be added to spectra)
    """

    id = models.CharField(primary_key=True, max_length=22, default=shortuuid, editable=False)
    extra_lights = models.ManyToManyField("Light", blank=True, related_name="microscopes")
    extra_cameras = models.ManyToManyField("Camera", blank=True, related_name="microscopes")
    extra_lasers = ArrayField(
        models.PositiveSmallIntegerField(validators=[MinValueValidator(300), MaxValueValidator(1600)]),
        default=list,
        blank=True,
    )
    collection = models.ForeignKey(
        "ProteinCollection",
        blank=True,
        null=True,
        related_name="on_scope",
        on_delete=models.CASCADE,
    )
    fluors = models.ForeignKey(
        "FluorophoreCollection",
        blank=True,
        null=True,
        related_name="fluor_on_scope",
        on_delete=models.CASCADE,
    )

    cfg_calc_efficiency = models.BooleanField(
        default=True,
        help_text="Calculate efficiency on update.",
    )
    cfg_fill_area = models.BooleanField(
        default=True,
        help_text="Fill area under spectra.",
    )
    cfg_min_wave = models.PositiveSmallIntegerField(
        default=350,
        validators=[MinValueValidator(300), MaxValueValidator(1199)],
        help_text="Minimum wavelength to display on page load.",
    )
    cfg_max_wave = models.PositiveSmallIntegerField(
        default=800,
        validators=[MinValueValidator(300), MaxValueValidator(1199)],
        help_text="Maximum wavelength to display on page load.",
    )
    cfg_enable_pan_zoom = models.BooleanField(
        default=True,
        help_text="Enable pan and zoom on spectra plot.",
    )

    class Meta:
        ordering = ["created"]

    @classmethod
    def from_oclist(cls, name, oclist):
        if not isinstance(oclist, list | tuple):
            raise ValueError("oclist argument must be list or tuple")
        microscope = cls(name=name)
        microscope.save()
        for row in oclist:
            newoc = quick_OC(row[0], row[1:], microscope)
            if newoc:
                microscope.optical_configs.add(newoc)
        microscope.save()
        return microscope

    @cached_property
    def has_inverted_bs(self):
        return self.optical_configs.filter(
            filterplacement__path=FilterPlacement.BS, filterplacement__reflects=True
        ).exists()

    @cached_property
    def has_reflective_emfilters(self):
        return self.optical_configs.filter(
            filterplacement__path=FilterPlacement.EM, filterplacement__reflects=True
        ).exists()

    def inverted_bs_set(self):
        return {i.filter.slug for oc in self.optical_configs.all() for i in oc.inverted_bs.all()}

    def get_absolute_url(self):
        return reverse("proteins:microscope-detail", args=[self.id])

    @cached_property
    def lights(self):
        oclights = Light.objects.filter(id__in=self.optical_configs.values("light"))
        return oclights

    @cached_property
    def cameras(self):
        occams = Camera.objects.filter(id__in=self.optical_configs.values("camera"))
        return occams

    @cached_property
    def lasers(self):
        return list(
            self.optical_configs.exclude(laser=None)
            .values_list("laser", flat=True)
            .order_by("laser")
            .distinct("laser")
        )

    @cached_property
    def ex_filters(self):
        return Filter.objects.filter(
            id__in=FilterPlacement.objects.filter(
                config__id__in=self.optical_configs.values("id"),
                path=FilterPlacement.EX,
            )
            .distinct("filter")
            .values("filter__id")
        )

    @cached_property
    def em_filters(self):
        return Filter.objects.filter(
            id__in=FilterPlacement.objects.filter(
                config__id__in=self.optical_configs.values("id"),
                path=FilterPlacement.EM,
            )
            .distinct("filter")
            .values("filter__id")
        )

    @cached_property
    def bs_filters(self):
        return Filter.objects.filter(
            id__in=FilterPlacement.objects.filter(
                config__id__in=self.optical_configs.values("id"),
                path=FilterPlacement.BS,
            )
            .distinct("filter")
            .values("filter__id")
        )

    @cached_property
    def spectra(self):
        spectra = []
        for f in (
            self.ex_filters,
            self.em_filters,
            self.bs_filters,
            self.lights,
            self.cameras,
        ):
            for i in f.select_related("spectrum"):
                spectra.append(i.spectrum)
        return spectra

    def spectra_d3(self):
        return [spec.d3dict() for spec in self.spectra]


def invert(sp):
    return [[a[0], 1 - a[1]] for a in sp]


OC_CACHE_KEY = "optical_config_list"


def get_cached_optical_configs(timeout=60 * 60):
    ocinfo = cache.get(OC_CACHE_KEY)
    if not ocinfo:
        vals = OpticalConfig.objects.all().values("id", "name", "comments", "microscope__id", "microscope__name")
        ocinfo = []
        for val in vals:
            scope = {
                "id": val.pop("microscope__id"),
                "name": val.pop("microscope__name"),
            }
            val["microscope"] = scope
            ocinfo.append(val)
        ocinfo = json.dumps({"data": {"opticalConfigs": ocinfo}})
        cache.set(OC_CACHE_KEY, ocinfo, timeout)
    return ocinfo


class OpticalConfig(OwnedCollection):
    """A a single optical configuration comprising a set of filters"""

    microscope = models.ForeignKey("Microscope", related_name="optical_configs", on_delete=models.CASCADE)
    comments = models.CharField(max_length=256, blank=True)
    filters = models.ManyToManyField("Filter", related_name="optical_configs", blank=True, through="FilterPlacement")
    light = models.ForeignKey(
        "Light",
        null=True,
        blank=True,
        related_name="optical_configs",
        on_delete=models.SET_NULL,
    )
    camera = models.ForeignKey(
        "Camera",
        null=True,
        blank=True,
        related_name="optical_configs",
        on_delete=models.SET_NULL,
    )
    laser = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(300), MaxValueValidator(1600)],
    )

    class Meta:
        unique_together = (("name", "microscope"),)
        ordering = ["name"]

    def save(self, **kwargs):
        cache.delete(OC_CACHE_KEY)
        super().save(**kwargs)

    @cached_property
    def ex_filters(self):
        """all filters that have an excitation role"""
        return self.filters.filter(filterplacement__path=FilterPlacement.EX)

    @cached_property
    def em_filters(self):
        """all filters that have an emission role"""
        return self.filters.filter(filterplacement__path=FilterPlacement.EM, filterplacement__reflects=False)

    @cached_property
    def bs_filters(self):
        """all filters that are in both ex and em paths have a beamsplitting role"""
        return self.filters.filter(filterplacement__path=FilterPlacement.BS)

    @cached_property
    def ref_em_filters(self):
        return self.filters.filter(filterplacement__path=FilterPlacement.EM, filterplacement__reflects=True)

    @cached_property
    def ex_spectra(self):
        """returns components in the excitation path"""
        p = []
        if self.laser:
            p.append([[self.laser - 1, 0], [self.laser, 1], [self.laser + 1, 0]])
        elif self.light:
            p.append(self.light.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.BS):
            p.append(invert(x.filter.spectrum.data) if not x.reflects else x.filter.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.EX):
            p.append(invert(x.filter.spectrum.data) if x.reflects else x.filter.spectrum.data)
        return p

    @cached_property
    def em_spectra(self):
        """returns components in the emissino path"""
        p = []
        if self.camera:
            p.append(self.camera.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.BS):
            p.append(invert(x.filter.spectrum.data) if x.reflects else x.filter.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.EM):
            p.append(invert(x.filter.spectrum.data) if x.reflects else x.filter.spectrum.data)
        return p

    @cached_property
    def combined_ex_spectra(self):
        return spectral_product(self.ex_spectra)

    @cached_property
    def combined_em_spectra(self):
        return spectral_product(self.em_spectra)

    @cached_property
    def inverted_bs(self):
        return self.filterplacement_set.filter(path=FilterPlacement.BS, reflects=True)

    def add_filter(self, filter, path, reflects=False):
        return FilterPlacement.objects.create(filter=filter, config=self, reflects=reflects, path=path)

    def __repr__(self):
        fltrs = sorted_ex2em(self.filters.all())
        return f"<{self.__class__.__name__}: {', '.join([f.name for f in fltrs])}>"

    def __str__(self):
        return super().__str__() or self.__repr__().lstrip("<").rstrip(">")

    def get_absolute_url(self):
        return "{}?c={}".format(
            reverse("proteins:microscope-detail", args=[self.microscope.id]),
            urllib.parse.quote(self.name),
        )


class FilterPlacement(models.Model):
    """Through table to specify placement of each filter in an OpticalConfig"""

    EX = "ex"
    EM = "em"
    BS = "bs"
    PATH_CHOICES = ((EX, "Excitation Path"), (EM, "Emission Path"), (BS, "Both Paths"))

    filter = models.ForeignKey("Filter", on_delete=models.CASCADE)
    config = models.ForeignKey("OpticalConfig", on_delete=models.CASCADE)
    path = models.CharField(max_length=2, choices=PATH_CHOICES, verbose_name="Ex/Bs/Em Path")
    # when path == BS, reflects refers to the emission path
    reflects = models.BooleanField(default=False, help_text="Filter reflects emission (if BS or EM filter)")

    def __str__(self):
        return self.__repr__().lstrip("<").rstrip(">")

    def __repr__(self):
        return f"<{self.path.title()} Filter: {self.filter.name}{' (reflecting)' if self.reflects else ''}>"


def quick_OC(name, filternames, scope, bs_ex_reflect=True):
    """generate an optical config from a tuple of strings (or tuples)

    valid examples:
        # tuple len 3
        quick_OC('Widefield Blue', ('ET395/25x', 'T425lpxr', 'ET460/50m'))
        # tuple len 4, the last value sets whether beamsplitter reflects exc.
        quick_OC('447 confocal', ('442', '445/515/561', 'ET480/40m', True))

        order of tuple matters:
            (excitation filters, beam splitters, emission filters)
        any integers in the first position will be assumed to be lasers
        tuples can be used to add multiples:
            quick_OC('double em', ['ET395/25x', 'T425lpxr', ('ET460/50m', 'ET605/70m')])
    """
    if len(filternames) == 4:
        bs_em_reflect = not bool(filternames.pop())
    if not len(filternames) == 3:
        raise ValueError("filternames argument must be iterable with length 3 or 4")

    oc = OpticalConfig.objects.create(name=name, microscope=scope)

    def _assign_filt(_f, i, iexact=False):
        # i = 0 is excitation
        # i = 1 is dichroic
        # i = 2 is emission
        _paths = [FilterPlacement.EX, FilterPlacement.BS, FilterPlacement.EM]
        try:
            if i != 0:
                raise ValueError()
            oc.laser = int(_f)
            oc.save()
        except ValueError:
            if iexact:
                filt = Filter.objects.get(part__iexact=_f)
            else:
                filt = Filter.objects.get(name__icontains=_f)
            fp = FilterPlacement(
                filter=filt,
                config=oc,
                path=_paths[i],
                reflects=bs_em_reflect if i == 1 else False,
            )
            fp.save()

    for i, fnames in enumerate(filternames):
        if not fnames:
            continue
        if isinstance(fnames, str):
            fnames = [fnames]
        elif isinstance(fnames, tuple | list):
            pass
        else:
            raise ValueError(f"value must be string, list, or tuple: {fnames}")

        for fname in fnames:
            try:
                _assign_filt(fname, i)
            except MultipleObjectsReturned:
                try:
                    _assign_filt(fname, i, True)
                except ObjectDoesNotExist:
                    print(f'Filter name "{fname}" returned multiple hits and exact match not found')
                    oc.delete()
                    return None
            except ObjectDoesNotExist:
                print(f'Filter name "{fname}" not found')
                oc.delete()
                return None
    return oc
