import graphene
import graphene_django_optimizer as gdo
from django.utils.text import slugify
from graphene_django.filter import DjangoFilterConnectionField
from graphql import FieldNode, GraphQLError, GraphQLResolveInfo

from proteins import models
from proteins.filters import ProteinFilter
from proteins.models.spectrum import get_spectra_list
from proteins.schema import relay, types


def get_spectrum(id):
    # (not cached here: the GraphQL view caches whole responses)
    return (
        models.Spectrum.objects.filter(id=id)
        .select_related("owner_fluor", "owner_camera", "owner_filter", "owner_light")
        .first()
    )


def get_requested_fields(info: GraphQLResolveInfo) -> set[str]:
    if not info.field_nodes or not (selection_set := info.field_nodes[0].selection_set):
        return set()
    return {node.name.value for node in selection_set.selections if isinstance(node, FieldNode)}


class Query(graphene.ObjectType):
    # this relay query delivers filterable paginated results
    all_proteins = DjangoFilterConnectionField(relay.ProteinNode, filterset_class=ProteinFilter)

    def resolve_all_proteins(self, info, **kwargs):
        # (a connection is not optimized automatically, unlike a list of proteins)
        return gdo.query(models.Protein.objects.all(), info)

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
            except models.Protein.DoesNotExist:
                return None
        return None

    # spectra = graphene.List(Spectrum)
    spectra = graphene.List(
        types.SpectrumInfo, subtype=graphene.String(), category=graphene.String()
    )
    spectrum = graphene.Field(types.Spectrum, id=graphene.Int())

    def resolve_spectra(self, info, **kwargs):
        requested_fields = get_requested_fields(info)

        fkwargs = {}
        if subtype := kwargs.get("subtype"):
            fkwargs["subtype"] = str(subtype).lower()
        if cat := kwargs.get("category"):
            fkwargs["category"] = str(cat).lower()

        if "owner" in requested_fields:
            # Use the optimized get_spectra_list function (no caching for GraphQL)
            return get_spectra_list(**fkwargs)
        elif fkwargs:
            return models.Spectrum.objects.filter(**fkwargs).values(*requested_fields)
        else:
            return models.Spectrum.objects.all().values(*requested_fields)

    def resolve_spectrum(self, info, **kwargs):
        _id = kwargs.get("id")
        return get_spectrum(_id) if _id is not None else None

    # def resolve_spectra(self, info, **kwargs):
    #     return gdo.query(models.Spectrum.objects.all(), info)

    state = graphene.Field(types.State, id=graphene.Int())
    states = graphene.List(types.State)

    def resolve_states(self, info, **kwargs):
        # (ordered: states have no ordering of their own, and an optimized query
        # returns them in another order than the plain one did)
        return gdo.query(models.State.objects.order_by("id"), info)

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

    dyes = graphene.List(types.DyeState)
    dye = graphene.Field(types.DyeState, id=graphene.Int(), name=graphene.String())

    def resolve_dyes(self, info, **kwargs):
        return gdo.query(models.DyeState.objects.all(), info)

    def resolve_dye(self, info, **kwargs):
        name = kwargs.get("name")
        if name is not None:
            slug = slugify(name)
            return gdo.query(models.DyeState.objects.filter(slug=slug), info).get()
        _id = kwargs.get("id")
        if _id is not None:
            return gdo.query(models.DyeState.objects.filter(id=_id), info).get()
        return None
