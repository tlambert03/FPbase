import graphene
from django.core.cache import cache
from graphene_django.filter import DjangoFilterConnectionField
from graphql import FieldNode, GraphQLError, GraphQLResolveInfo

from .. import models
from ..filters import ProteinFilter
from . import _optimizer as gdo
from . import relay, types


def get_cached_spectrum(id, timeout=60 * 60 * 24):
    key = f"_spectrum_{id}"
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


def get_requested_fields(info: GraphQLResolveInfo) -> set[str]:
    if not info.field_nodes or not (selection_set := info.field_nodes[0].selection_set):
        return set()
    return {node.name.value for node in selection_set.selections if isinstance(node, FieldNode)}


class Query(graphene.ObjectType):
    # this relay query delivers filterable paginated results
    all_proteins = DjangoFilterConnectionField(relay.ProteinNode, filterset_class=ProteinFilter)

    microscopes = graphene.List(types.Microscope)
    microscope = graphene.Field(types.Microscope, id=graphene.String())

    def resolve_microscopes(self, info, **kwargs):
        return gdo.query(models.Microscope.objects.all(), info)

    def resolve_microscope(self, info, **kwargs):
        _id = kwargs.get("id")
        if _id is not None:
            try:
                obj = gdo.query(models.Microscope.objects.filter(id__istartswith=_id), info)
                return obj.get()
            except models.Microscope.MultipleObjectsReturned as e:
                raise GraphQLError(f'Multiple microscopes found starting with "{_id}"') from e
            except models.Microscope.DoesNotExist:
                return None
        return None

    organisms = graphene.List(types.Organism)
    organism = graphene.Field(types.Organism, id=graphene.Int())

    def resolve_organisms(self, info, **kwargs):
        return gdo.query(models.Organism.objects.all(), info)

    def resolve_organism(self, info, **kwargs):
        _id = kwargs.get("id")
        if _id is not None:
            try:
                return gdo.query(models.Organism.objects.filter(id=_id), info).get()
            except models.Organism.DoesNotExist:
                return None
        return None

    proteins = graphene.List(types.Protein)
    protein = graphene.Field(types.Protein, id=graphene.String())

    def resolve_proteins(self, info, **kwargs):
        return gdo.query(models.Protein.objects.all(), info)

    def resolve_protein(self, info, **kwargs):
        _id = kwargs.get("id")
        if _id is not None:
            try:
                return gdo.query(models.Protein.objects.filter(uuid=_id), info).get()
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
        _id = kwargs.get("id")
        return get_cached_spectrum(_id) if _id is not None else None

    # def resolve_spectra(self, info, **kwargs):
    #     return gdo.query(models.Spectrum.objects.all(), info)

    state = graphene.Field(types.State, id=graphene.Int())
    states = graphene.List(types.State)

    def resolve_states(self, info, **kwargs):
        return models.State.objects.all()

    def resolve_state(self, info, **kwargs):
        _id = kwargs.get("id")
        return models.State.objects.get(id=_id) if _id is not None else None

    opticalConfigs = graphene.List(types.OpticalConfig)
    opticalConfig = graphene.Field(types.OpticalConfig, id=graphene.Int())

    def resolve_opticalConfigs(self, info, **kwargs):
        # return models.OpticalConfig.objects.all().prefetch_related("microscope")
        return gdo.query(models.OpticalConfig.objects.all(), info)

    def resolve_opticalConfig(self, info, **kwargs):
        _id = kwargs.get("id")
        if _id is not None:
            return gdo.query(models.OpticalConfig.objects.filter(id=_id), info).get()
        return None
