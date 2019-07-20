import graphene
import graphene_django_optimizer as gdo
from graphql import GraphQLError
from graphene_django.filter import DjangoFilterConnectionField

from . import types, relay
from .. import models
from ..filters import ProteinFilter

from django.core.cache import cache


def get_cached_spectrum(id, timeout=60 * 60 * 24):
    key = "_spectrum_{}".format(id)
    spectrum = cache.get(key)
    if not spectrum:
        try:
            spectrum = (
                models.Spectrum.objects.filter(id=id)
                .select_related(
                    "owner_state",
                    "owner_state__protein",
                    "owner_dye",
                    "owner_camera",
                    "owner_filter",
                    "owner_light",
                )
                .get()
            )
            cache.set(key, spectrum, timeout)
        except models.Spectrum.DoesNotExist:
            return None
    return spectrum


def get_requested_fields(info):
    selections = info.field_asts[0].selection_set.selections
    requested_fields = [f.name.value for f in selections]
    return requested_fields


class Query(graphene.ObjectType):

    # this relay query delivers filterable paginated results
    all_proteins = DjangoFilterConnectionField(
        relay.ProteinNode, filterset_class=ProteinFilter
    )

    microscopes = graphene.List(types.Microscope)
    microscope = graphene.Field(types.Microscope, id=graphene.String())

    def resolve_microscopes(self, info, **kwargs):
        return gdo.query(models.Microscope.objects.all(), info)

    def resolve_microscope(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            try:
                obj = gdo.query(
                    models.Microscope.objects.filter(id__istartswith=id), info
                )
                return obj.get()
            except models.Microscope.MultipleObjectsReturned:
                raise GraphQLError(
                    'Multiple microscopes found starting with "{}"'.format(id)
                )
            except models.Microscope.DoesNotExist:
                return None
        return None

    organisms = graphene.List(types.Organism)
    organism = graphene.Field(types.Organism, id=graphene.Int())

    def resolve_organisms(self, info, **kwargs):
        return gdo.query(models.Organism.objects.all(), info)

    def resolve_organism(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            try:
                return gdo.query(models.Organism.objects.filter(id=id), info).get()
            except models.Organism.DoesNotExist:
                return None
        return None

    proteins = graphene.List(types.Protein)
    protein = graphene.Field(types.Protein, id=graphene.String())

    def resolve_proteins(self, info, **kwargs):
        return gdo.query(models.Protein.objects.all(), info)

    def resolve_protein(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            try:
                return gdo.query(models.Protein.objects.filter(uuid=id), info).get()
            except models.Spectrum.DoesNotExist:
                return None
        return None

    # spectra = graphene.List(Spectrum)
    spectra = graphene.List(types.SpectrumInfo)
    spectrum = graphene.Field(types.Spectrum, id=graphene.Int())

    def resolve_spectra(self, info, **kwargs):
        requested_fields = get_requested_fields(info)
        if "owner" in requested_fields:
            return models.Spectrum.objects.sluglist()
        return models.Spectrum.objects.all().values(*requested_fields)

    def resolve_spectrum(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return get_cached_spectrum(id)
        return None

    # def resolve_spectra(self, info, **kwargs):
    #     return gdo.query(models.Spectrum.objects.all(), info)

    state = graphene.Field(types.State, id=graphene.Int())
    states = graphene.List(types.State)

    def resolve_states(self, info, **kwargs):
        return models.State.objects.all()

    def resolve_state(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return models.State.objects.get(id=id)
        return None

    opticalConfigs = graphene.List(types.OpticalConfig)
    opticalConfig = graphene.Field(types.OpticalConfig, id=graphene.Int())

    def resolve_opticalConfigs(self, info, **kwargs):
        # return models.OpticalConfig.objects.all().prefetch_related("microscope")
        return gdo.query(models.OpticalConfig.objects.all(), info)

    def resolve_opticalConfig(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return gdo.query(models.OpticalConfig.objects.filter(id=id), info).get()
        return None
