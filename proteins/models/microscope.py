from django.db import models
from django.db.models import Count, Q
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse
from django.core.validators import MaxValueValidator, MinValueValidator
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from .spectrum import Filter
from .spectrum import sorted_ex2em
from ..util.helpers import shortuuid
from .collection import OwnedCollection


class Microscope(OwnedCollection):
    """ A microscope or other collection of optical configurations

    Items are intended to be added as configs, but can be directly added as
    filters/cameras/lights etc as well.

        configs: optical configs in the collection
                 (all spectra will be added to spectra)
    """
    id = models.CharField(primary_key=True, max_length=22, default=shortuuid, editable=False)
    ex_filters = models.ManyToManyField('Filter', blank=True, related_name='as_ex_filter')
    bs_filters = models.ManyToManyField('Filter', blank=True, related_name='as_bs_filter')
    em_filters = models.ManyToManyField('Filter', blank=True, related_name='as_em_filter')
    lights = models.ManyToManyField('Light', blank=True, related_name='microscopes')
    cameras = models.ManyToManyField('Camera', blank=True, related_name='microscopes')
    lasers = ArrayField(models.PositiveSmallIntegerField(
                        validators=[MinValueValidator(300), MaxValueValidator(1600)]),
                        default=list, blank=True)
    collection = models.ForeignKey('ProteinCollection', blank=True, null=True, related_name='on_scope', on_delete=models.CASCADE)
    fluors = models.ForeignKey('FluorophoreCollection', blank=True, null=True, related_name='fluor_on_scope', on_delete=models.CASCADE)

    @classmethod
    def from_oclist(cls, name, oclist):
        if not isinstance(oclist, (list, tuple)):
            raise ValueError('oclist argument must be list or tuple')
        microscope = cls(name=name)
        microscope.save()
        for row in oclist:
            newoc = quick_OC(row[0], row[1:], microscope)
            if newoc:
                microscope.optical_configs.add(newoc)
        microscope.save()
        return microscope

    def has_inverted_bs(self):
        for oc in self.optical_configs.all():
            if oc.inverted_bs.exists():
                return True
        return False

    def inverted_bs_set(self):
        return {i.filter.slug for oc in self.optical_configs.all() for i in oc.inverted_bs.all()}

    def save(self, *args):
        # add any filterset members to the spectra set
        if self.pk and self.optical_configs:
            # For now... this makes it so that a microscope can only have filters
            # that are used in at least one optical config
            self.ex_filters.clear()
            self.em_filters.clear()
            self.bs_filters.clear()
            # self.lights.clear()
            # self.cameras.clear()
            # self.lasers = []
            for oc in self.optical_configs.all():
                for _fl in ('ex_filters', 'bs_filters', 'em_filters'):
                    [getattr(self, _fl).add(x) for x in getattr(oc, _fl).all()]
                if oc.light:
                    self.lights.add(oc.light)
                if oc.camera:
                    self.cameras.add(oc.camera)
                if oc.laser:
                    self.lasers.append(oc.laser)
        self.lasers = list(set(self.lasers))
        self.full_clean()
        return super().save(*args)

    def get_absolute_url(self):
        return reverse("proteins:microscope-detail", args=[self.id])

    @property
    def spectra(self):
        spectra = []
        for f in (self.ex_filters, self.em_filters, self.bs_filters,
                  self.lights, self.cameras):
            for i in f.select_related('spectrum'):
                spectra.append(i.spectrum)
        return spectra

    def spectra_d3(self):
        return [spec.d3dict() for spec in self.spectra]


def invert(sp):
    return [[a[0], 1 - a[1]] for a in sp]


class OpticalConfig(OwnedCollection):
    """ A a single optical configuration comprising a set of filters """
    microscope = models.ForeignKey('Microscope', related_name='optical_configs', on_delete=models.CASCADE)
    filters = models.ManyToManyField(
        'Filter', related_name='optical_configs', blank=True, through='FilterPlacement')
    light = models.ForeignKey('Light', null=True, blank=True, related_name='optical_configs', on_delete=models.CASCADE)
    camera = models.ForeignKey('Camera', null=True, blank=True, related_name='optical_configs', on_delete=models.CASCADE)
    laser = models.PositiveSmallIntegerField(
        blank=True, null=True,
        validators=[MinValueValidator(300), MaxValueValidator(1600)])

    class Meta:
        unique_together = (("name", "microscope"),)
        ordering = ['name']

    def save(self, **kwargs):
        super().save(**kwargs)

    @property
    def ex_filters(self):
        """ all filters that have an excitation role """
        x = self.filters.filter(filterplacement__path=FilterPlacement.EX)
        return x.exclude(id__in=self.bs_filters)

    @property
    def em_filters(self):
        """ all filters that have an emission role """
        m = self.filters.filter(filterplacement__path=FilterPlacement.EM)
        return m.exclude(id__in=self.bs_filters)

    @property
    def bs_filters(self):
        """ all filters that are in both ex and em paths have a beamsplitting role """
        inx = Count('filterplacement__path', distinct=True,
                    filter=Q(filterplacement__path=FilterPlacement.EX))
        inm = Count('filterplacement__path',
                    filter=Q(filterplacement__path=FilterPlacement.EM))
        return self.filters.annotate(nx=inx, nm=inm).filter(nx__gt=0, nm__gt=0)

    @property
    def ex_spectra(self):
        """ returns components in the excitation path """
        p = []
        if self.laser:
            p.append([[self.laser - 1, 0], [self.laser, 1], [self.laser + 1, 0]])
        elif self.light:
            p.append(self.light.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.EX):
            p.append(invert(x.filter.spectrum.data) if x.reflects
                     else x.filter.spectrum.data)
        return p

    @property
    def em_spectra(self):
        """ returns components in the emissino path """
        p = []
        if self.camera:
            p.append(self.camera.spectrum.data)
        for x in self.filterplacement_set.filter(path=FilterPlacement.EM):
            p.append(invert(x.filter.spectrum.data) if x.reflects
                     else x.filter.spectrum.data)
        return p

    @property
    def inverted_bs(self):
        return self.filterplacement_set.filter(path=FilterPlacement.EM, reflects=True)

    def add_em_filter(self, filter, reflects=False):
        fp = FilterPlacement(filter=filter, config=self, reflects=reflects,
                             path=FilterPlacement.EM)
        fp.save()

    def add_ex_filter(self, filter, reflects=False):
        fp = FilterPlacement(filter=filter, config=self, reflects=reflects,
                             path=FilterPlacement.EX)
        fp.save()

    def add_bs_filter(self, filter, em_reflects=False):
        fpx = FilterPlacement(filter=filter, config=self,
                              reflects=not em_reflects,
                              path=FilterPlacement.EX)
        fpm = FilterPlacement(filter=filter, config=self,
                              reflects=em_reflects,
                              path=FilterPlacement.EM)
        fpx.save()
        fpm.save()

    def __repr__(self):
        fltrs = sorted_ex2em(self.filters.all())
        return "<{}: {}>".format(
            self.__class__.__name__, ", ".join([f.name for f in fltrs]))

    def __str__(self):
        return super().__str__() or self.__repr__().lstrip('<').rstrip('>')

    def get_absolute_url(self):
        return reverse("proteins:filterset-detail", args=[self.id])


class FilterPlacement(models.Model):
    """ Through table to specify placement of each filter in an OpticalConfig """

    EX = 'ex'
    EM = 'em'
    PATH_CHOICES = (
        (EX, 'Excitation Path'),
        (EM, 'Emission Path'),
    )

    filter = models.ForeignKey('Filter', on_delete=models.CASCADE)
    config = models.ForeignKey('OpticalConfig', on_delete=models.CASCADE)
    path = models.CharField(max_length=2, choices=PATH_CHOICES, verbose_name='Ex/Em Path')
    reflects = models.BooleanField(
        default=False,
        help_text='Filter reflects light at this position in the light path')

    def __repr__(self):
        return "<{} Filter: {}{}>".format(self.path.title(), self.filter.name,
                                          " (reflecting)" if self.reflects else '')


def quick_OC(name, filternames, scope, bs_ex_reflect=True):
    """ generate an optical config from a tuple of strings (or tuples)

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
        bs_ex_reflect = bool(filternames.pop())
    if not len(filternames) == 3:
        raise ValueError('filternames argument must be iterable with length 3 or 4')

    oc = OpticalConfig.objects.create(name=name, microscope=scope)

    def _assign_filt(_f, i, iexact=False):
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
            if i == 1:
                fpx = FilterPlacement(filter=filt, config=oc,
                                      reflects=bs_ex_reflect,
                                      path=FilterPlacement.EX)
                fpm = FilterPlacement(filter=filt, config=oc,
                                      reflects=not bs_ex_reflect,
                                      path=FilterPlacement.EM)
                fpx.save()
                fpm.save()
            else:
                _path = FilterPlacement.EX if i < 1 else FilterPlacement.EM
                fp = FilterPlacement(filter=filt, config=oc, path=_path)
                fp.save()

    for i, fnames in enumerate(filternames):
        if not fnames:
            continue
        if isinstance(fnames, str):
            fnames = [fnames]
        elif isinstance(fnames, (tuple, list)):
            pass
        else:
            raise ValueError('value must be string, list, or tuple: %s' % fnames)

        for fname in fnames:
            try:
                _assign_filt(fname, i)
            except MultipleObjectsReturned:
                try:
                    _assign_filt(fname, i, True)
                except ObjectDoesNotExist:
                    print('Filter name "%s" returned multiple hits and exact match not found' % fname)
                    oc.delete()
                    return None
            except ObjectDoesNotExist:
                print('Filter name "%s" not found' % fname)
                oc.delete()
                return None
    return oc
