import graphene
import graphene_django_optimizer as gdo
from django.db.models import Prefetch
from graphene.utils.str_converters import to_camel_case
from graphene_django.converter import get_choices
from graphene_django.types import DjangoObjectType
from graphql import GraphQLError
from graphene_django.filter import DjangoFilterConnectionField


from references.schema import Reference

from . import models


def nullable_enum_from_field(_model, _field):
    field = _model._meta.get_field(_field)
    choices = getattr(field, "choices", None)
    if choices:
        meta = field.model._meta
        name = to_camel_case("my{}_{}".format(meta.object_name, field.name))
        choices = list(get_choices(choices))
        named_choices = [(c[0], c[1]) for c in choices]
        named_choices_descriptions = {c[0]: c[2] for c in choices}

        class EnumWithDescriptionsType(object):
            @property
            def description(self):
                return named_choices_descriptions[self.name]

        enum = graphene.Enum(name, list(named_choices), type=EnumWithDescriptionsType)
        converted = enum(
            description=field.help_text, required=not (field.null or field.blank)
        )
    else:
        raise NotImplementedError("Field does NOT have choices")
    return converted


def parse_selection(sel):
    return (
        {sel.name.value: [parse_selection(s) for s in sel.selection_set.selections]}
        if sel.selection_set
        else sel.name.value
    )


def get_requested_fields(info):
    selections = info.field_asts[0].selection_set.selections
    requested_fields = [f.name.value for f in selections]
    return requested_fields


class Organism(gdo.OptimizedDjangoObjectType):
    proteins = graphene.List(lambda: Protein)

    class Meta:
        model = models.Organism

    @gdo.resolver_hints(select_related=("proteins"), only=("proteins"))
    def resolve_proteins(self, info):
        return self.proteins.all()


class OSERMeasurement(gdo.OptimizedDjangoObjectType):
    class Meta:
        model = models.OSERMeasurement


class StateTransition(gdo.OptimizedDjangoObjectType):
    fromState = graphene.Field(lambda: State)
    toState = graphene.Field(lambda: State)

    class Meta:
        model = models.StateTransition

    @gdo.resolver_hints(select_related=("from_state"), only=("from_state"))
    def resolve_fromState(self, info):
        return self.from_state

    @gdo.resolver_hints(select_related=("to_state"), only=("to_state"))
    def resolve_toState(self, info):
        return self.to_state


class Protein(gdo.OptimizedDjangoObjectType):
    id = graphene.ID(source="uuid", required=True)
    parentOrganism = graphene.Field(Organism)
    primaryReference = graphene.Field(Reference)
    switchType = nullable_enum_from_field(models.Protein, "switch_type")
    agg = nullable_enum_from_field(models.Protein, "agg")
    cofactor = nullable_enum_from_field(models.Protein, "cofactor")
    oser = graphene.List(OSERMeasurement)
    transitions = graphene.List(StateTransition)

    class Meta:
        model = models.Protein
        exclude = (
            "id",
            "status",
            "status_changed",
            "uuid",
            "base_name",
            "switch_type",
            "agg",
            "cofactor",
            "oser",
            "transitions",
        )

    @gdo.resolver_hints(
        prefetch_related=lambda info: Prefetch(
            "transitions",
            queryset=gdo.query(models.StateTransition.objects.all(), info),
        )
    )
    def resolve_transitions(self, info, **kwargs):
        return self.transitions.all()

    def resolve_switchType(self, info):
        return self.switch_type or None

    @gdo.resolver_hints(prefetch_related=("oser_measurements"))
    def resolve_oser(self, info):
        return self.oser_measurements.all()

    def resolve_agg(self, info):
        return self.agg or None

    def resolve_cofactor(self, info):
        return self.cofactor or None

    @gdo.resolver_hints(select_related=("parent_organism"), only=("parent_organism"))
    def resolve_parentOrganism(self, info):
        return self.parent_organism

    @gdo.resolver_hints(
        select_related=("primary_reference"), only=("primary_reference")
    )
    def resolve_primaryReference(self, info):
        return self.primary_reference


class ProteinNode(Protein):
    class Meta:
        model = models.Protein
        interfaces = (graphene.relay.Node,)
        filter_fields = {
            "name": ["exact", "icontains", "istartswith"],
            "switch_type": ["exact"],
            "agg": ["exact"],
        }


class FluorophoreInterface(graphene.Interface):
    qy = graphene.Float()
    extCoeff = graphene.Float(source="ext_coeff")
    twopPeakgm = graphene.Float(source="twop_peakGM")
    exMax = graphene.Float(source="ex_max")
    emMax = graphene.Float(source="em_max")


class SpectrumOwnerInterface(graphene.Interface):
    name = graphene.String()
    id = graphene.ID()
    typ = graphene.String()
    slug = graphene.String()

    def resolve_typ(self, info):
        return self.__class__.__name__.lower()

    def resolve_name(self, info):
        return str(self)


class Camera(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Camera


class Dye(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface, FluorophoreInterface)
        model = models.Dye


class Filter(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Filter


class Light(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Light


class State(gdo.OptimizedDjangoObjectType):
    protein = graphene.Field(Protein)

    class Meta:
        interfaces = (SpectrumOwnerInterface, FluorophoreInterface)
        model = models.State

    @gdo.resolver_hints(select_related=("protein",), only=("protein",))
    def resolve_protein(self, info, **kwargs):
        return self.protein

    # spectra = graphene.List(SpectrumType)

    # def resolve_spectra(self, info, **kwargs):
    #     return self.spectrumowner.spectra.all()


class SpectrumOwnerUnion(graphene.Union):
    class Meta:
        types = (State,)


class Spectrum(gdo.OptimizedDjangoObjectType):
    class Meta:
        model = models.Spectrum

    owner = graphene.Field(SpectrumOwnerInterface)
    color = graphene.String()

    @gdo.resolver_hints(
        select_related=(
            "owner_state",
            "owner_dye",
            "owner_camera",
            "owner_filter",
            "owner_light",
        ),
        only=(
            "owner_state",
            "owner_dye",
            "owner_camera",
            "owner_filter",
            "owner_light",
        ),
    )
    def resolve_owner(self, info, **kwargs):
        return self.owner

    def resolve_color(self, info, **kwargs):
        return self.color()


class SpectrumOwnerInfo(graphene.ObjectType):
    name = graphene.String()
    slug = graphene.String()
    url = graphene.String()
    id = graphene.ID()

    def resolve_name(self, info):
        return self.get("name")

    def resolve_slug(self, info):
        return self.get("slug")

    def resolve_url(self, info):
        return self.get("url")

    def resolve_id(self, info):
        return int(self.get("id"))


class SpectrumInfo(graphene.ObjectType):
    id = graphene.ID()
    category = graphene.String()
    subtype = graphene.String()
    owner = graphene.Field(SpectrumOwnerInfo)

    def resolve_id(self, info):
        return int(self.get("id"))

    def resolve_category(self, info):
        return self.get("category").upper()

    def resolve_subtype(self, info):
        return self.get("subtype").upper()

    def resolve_owner(self, info):
        return self.get("owner")


class FilterPlacement(gdo.OptimizedDjangoObjectType):
    id = graphene.ID()
    spectrum = graphene.Field(Spectrum)
    spectrumId = id = graphene.ID()
    reflects = graphene.Boolean()

    class Meta:
        model = models.FilterPlacement

    @gdo.resolver_hints(select_related=("filter__spectrum"))
    def resolve_spectrum(self, info):
        return self.filter.spectrum

    @gdo.resolver_hints(
        select_related=("filter__spectrum"), only=("filter__spectrum__id",)
    )
    def resolve_spectrumId(self, info):
        return self.filter.spectrum.id

    def resolve_id(self, info):
        return self.filter_id


class Microscope(DjangoObjectType):
    class Meta:
        model = models.Microscope
        only_fields = (
            "id",
            "name",
            "description",
            "extra_lights",
            "extra_cameras",
            "extra_lasers",
            "optical_configs",
        )


class OpticalConfig(gdo.OptimizedDjangoObjectType):
    filters = graphene.List(FilterPlacement)
    microscope = graphene.Field(Microscope)

    class Meta:
        model = models.OpticalConfig
        only_fields = (
            "id",
            "name",
            "description",
            "microscope",
            "comments",
            "filters",
            "light",
            "camera",
            "laser",
        )

    @gdo.resolver_hints(
        prefetch_related=(
            "filterplacement_set",
            "filterplacement_set__filter",
            "filterplacement_set__filter__spectrum",
        )
    )
    def resolve_filters(self, info):
        return self.filterplacement_set.all()


class Query(graphene.ObjectType):

    microscopes = graphene.List(Microscope)
    microscope = graphene.Field(Microscope, id=graphene.String())

    all_proteins = DjangoFilterConnectionField(ProteinNode)

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

    organisms = graphene.List(Organism)
    organism = graphene.Field(Organism, id=graphene.Int())

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

    proteins = graphene.List(Protein)
    protein = graphene.Field(Protein, id=graphene.String())

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
    spectra = graphene.List(SpectrumInfo)
    spectrum = graphene.Field(Spectrum, id=graphene.Int())

    def resolve_spectra(self, info, **kwargs):
        requested_fields = get_requested_fields(info)
        if "owner" in requested_fields:
            return models.Spectrum.objects.sluglist()
        return models.Spectrum.objects.all().values(*requested_fields)

    def resolve_spectrum(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            try:
                return gdo.query(models.Spectrum.objects.filter(id=id), info).get()
            except models.Spectrum.DoesNotExist:
                return None
        return None

    state = graphene.Field(State, id=graphene.Int())
    states = graphene.List(State)

    def resolve_states(self, info, **kwargs):
        return models.State.objects.all()

    def resolve_state(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return models.State.objects.get(id=id)
        return None

    opticalConfigs = graphene.List(OpticalConfig)
    opticalConfig = graphene.Field(OpticalConfig, id=graphene.Int())

    def resolve_opticalConfigs(self, info, **kwargs):
        # return models.OpticalConfig.objects.all().prefetch_related("microscope")
        return gdo.query(models.OpticalConfig.objects.all(), info)

    def resolve_opticalConfig(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return gdo.query(models.OpticalConfig.objects.filter(id=id), info).get()
        return None

    # def resolve_spectra(self, info, **kwargs):
    #     return gdo.query(models.Spectrum.objects.all(), info)
