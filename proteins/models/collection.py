from django.db import models
from django.db.models import Count, Q
from django.contrib.auth import get_user_model
from django.urls import reverse
from django.core.validators import MaxValueValidator, MinValueValidator
from model_utils.models import TimeStampedModel
from .spectrum import sorted_ex2em
from collections import namedtuple


User = get_user_model()

oriented_filter = namedtuple('OrientedFilter', ['filter', 'reflects'])


class OwnedCollection(TimeStampedModel):
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=512, blank=True)
    owner = models.ForeignKey(User, blank=True, null=True,
                              related_name='%(class)s' + 's',
                              on_delete=models.SET_NULL,)
    private = models.BooleanField(default=False, verbose_name="Private Collection",
                                  help_text="Private collections can not be seen "
                                            "by or shared with other users")

    def __str__(self):
        return self.name

    class Meta:
        abstract = True
        unique_together = (('owner', 'name'),)


class ProteinCollection(OwnedCollection):
    proteins = models.ManyToManyField('Protein', related_name='collection_memberships')

    def get_absolute_url(self):
        return reverse("proteins:collection-detail", args=[self.id])


class SpectrumCollection(OwnedCollection):
    """ A microscope or other collection of optical configurations

        spectra: contains ALL spectra in the collection
        configs: optical configs in the collection
                 (all spectra will be added to spectra)
    """
    spectra = models.ManyToManyField('Spectrum', related_name='collection_memberships')
    configs = models.ManyToManyField('OpticalConfig', related_name='collection_memberships')

    def __repr__(self):
        count = self.spectra.all().count()
        return "<{}: with {} spectra>".format(self.__class__.__name__, count)

    def save(self, *args):
        # add any filterset members to the spectra set
        if self.pk and self.configs:
            for oc in self.configs.all():
                for spec in oc.spectra_set:
                    self.spectra.add(spec.spectrum)
        return super().save(*args)

    def get_absolute_url(self):
        return reverse("proteins:spectrum-collection-detail", args=[self.id])


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


class OpticalConfig(OwnedCollection):
    """ A a single optical configuration comprising a set of filters """

    filters = models.ManyToManyField(
        'Filter', related_name='optical_configs', through='FilterPlacement')
    light = models.ForeignKey('Light', null=True, related_name='optical_configs', on_delete=models.CASCADE)
    camera = models.ForeignKey('Camera', null=True, related_name='optical_configs', on_delete=models.CASCADE)
    laser = models.PositiveSmallIntegerField(
        blank=True, null=True,
        validators=[MinValueValidator(300), MaxValueValidator(1600)])

    @property
    def ex_path(self):
        """ returns components in the excitation path """
        p = []
        if self.light:
            p.append((self.light, False))
        if self.laser:
            p.append((self.laser, False))
        for x in self.filterplacement_set.filter(path=FilterPlacement.EX):
            p.append(oriented_filter(x.filter, x.reflects))
        return p

    @property
    def em_path(self):
        """ returns components in the emissino path """
        p = []
        if self.camera:
            p.append((self.camera, False))
        for x in self.filterplacement_set.filter(path=FilterPlacement.EM):
            p.append(oriented_filter(x.filter, x.reflects))
        return p

    @property
    def both_paths(self):
        inx = Count('filterplacement__path',
                    filter=Q(filterplacement__path=FilterPlacement.EX))
        inm = Count('filterplacement__path',
                    filter=Q(filterplacement__path=FilterPlacement.EM))
        return self.filters.annotate(nx=inx, nm=inm).filter(nx__gt=0, nm__gt=0)

    @property
    def filter_set(self):
        """ returns unique set of all filters in this optical config """
        return set(self.filters.all())

    @property
    def spectra_set(self):
        return self.filter_set.union({s for s in (self.camera, self.light) if s})

    def __repr__(self):
        fltrs = sorted_ex2em(self.filters.all())
        return "<{}: {}>".format(
            self.__class__.__name__, ", ".join([f.name for f in fltrs]))

    def __str__(self):
        return super().__str__() or self.__repr__().lstrip('<').rstrip('>')

    def get_absolute_url(self):
        return reverse("proteins:filterset-detail", args=[self.id])


def quick_OC(name, filternames, bs_ex_reflect=True, bs_pos='both'):
    from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
    from proteins.models import Filter
    oc = OpticalConfig.objects.create(name=name)
    for i, fname in enumerate(filternames):
        if fname.isdigit():
            oc.laser = int(fname)
            oc.save()
            continue
        try:
            filt = Filter.objects.get(name__icontains=fname)
        except MultipleObjectsReturned:
            print('Filter name "%s" returned multiple hits' % fname)
            oc.delete()
            return None
        except ObjectDoesNotExist:
            print('Filter name "%s" not found' % fname)
            oc.delete()
            return None
        # put dichroics in both paths by default
        if filt and filt.subtype in (Filter.BS, Filter.LP):
            if bs_pos in ('ex', 'both'):
                fp = FilterPlacement(filter=filt, config=oc,
                                     reflects=bs_ex_reflect,
                                     path=FilterPlacement.EX)
                fp.save()
            if bs_pos in ('em', 'both'):
                fp = FilterPlacement(filter=filt, config=oc,
                                     reflects=not bs_ex_reflect,
                                     path=FilterPlacement.EM)
                fp.save()
        else:
            _path = FilterPlacement.EM if i else FilterPlacement.EX
            fp = FilterPlacement(filter=filt, config=oc, path=_path)
            fp.save()
    return oc
