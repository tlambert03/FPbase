import graphene
from graphene.utils.str_converters import to_camel_case
from graphene_django.converter import get_choices
from graphene_django.types import DjangoObjectType as BaseDjangoObjectType
from query_optimizer import DjangoListField, DjangoObjectType, RelatedField

from references.schema import Reference

from .. import models


def nullable_enum_from_field(_model, _field):
    field = _model._meta.get_field(_field)
    choices = getattr(field, "choices", None)
    if choices:
        meta = field.model._meta
        name = to_camel_case(f"my{meta.object_name}_{field.name}")
        choices = list(get_choices(choices))
        named_choices = [(c[0], c[1]) for c in choices]
        named_choices_descriptions = {c[0]: c[2] for c in choices}

        class EnumWithDescriptionsType:
            @property
            def description(self):
                return named_choices_descriptions[self.name]

        enum = graphene.Enum(name, list(named_choices), type=EnumWithDescriptionsType)
        converted = enum(description=field.help_text, required=not (field.null or field.blank))
    else:
        raise NotImplementedError("Field does NOT have choices")
    return converted


def parse_selection(sel):
    return (
        {sel.name.value: [parse_selection(s) for s in sel.selection_set.selections]}
        if sel.selection_set
        else sel.name.value
    )


class Organism(DjangoObjectType):
    proteins = DjangoListField(lambda: Protein)

    class Meta:
        model = models.Organism
        fields = "__all__"


class OSERMeasurement(DjangoObjectType):
    class Meta:
        model = models.OSERMeasurement
        fields = "__all__"


class StateTransition(DjangoObjectType):
    fromState = RelatedField(lambda: State, field_name="from_state")
    toState = RelatedField(lambda: State, field_name="to_state")

    class Meta:
        model = models.StateTransition
        fields = "__all__"


class Protein(DjangoObjectType):
    id = graphene.ID(source="uuid", required=True)
    parentOrganism = RelatedField(Organism, field_name="parent_organism")
    primaryReference = RelatedField(Reference, field_name="primary_reference")
    switchType = nullable_enum_from_field(models.Protein, "switch_type")
    agg = nullable_enum_from_field(models.Protein, "agg")
    cofactor = nullable_enum_from_field(models.Protein, "cofactor")
    oser = DjangoListField(OSERMeasurement, field_name="oser_measurements")
    transitions = DjangoListField(StateTransition)

    class Meta:
        model = models.Protein
        exclude = ("status", "status_changed", "uuid", "base_name", "switch_type")

    def resolve_switchType(self, info):
        return self.switch_type or None

    def resolve_agg(self, info):
        return self.agg or None

    def resolve_cofactor(self, info):
        return self.cofactor or None


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


class Camera(BaseDjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Camera
        fields = "__all__"


class Dye(BaseDjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface, FluorophoreInterface)
        model = models.Dye
        fields = "__all__"


class Filter(BaseDjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Filter
        fields = "__all__"


class Light(BaseDjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Light
        fields = "__all__"


class State(DjangoObjectType):
    protein = RelatedField(Protein)

    class Meta:
        interfaces = (SpectrumOwnerInterface, FluorophoreInterface)
        model = models.State
        fields = "__all__"

    # spectra = graphene.List(SpectrumType)

    # def resolve_spectra(self, info, **kwargs):
    #     return self.spectrumowner.spectra.all()


class SpectrumOwnerUnion(graphene.Union):
    class Meta:
        types = (State,)


# class SpectrumLoader(DataLoader):
#     def batch_load_fn(self, keys):
#         spectra = {
#             spectrum.id: spectrum
#             for spectrum in models.Spectrum.objects.filter(id__in=keys).select_related(
#                 "owner_state",
#                 "owner_state__protein",
#                 "owner_dye",
#                 "owner_camera",
#                 "owner_filter",
#                 "owner_light",
#             )
#         }
#         return Promise.resolve([spectra.get(spectra_id) for spectra_id in keys])


class Spectrum(DjangoObjectType):
    class Meta:
        model = models.Spectrum
        fields = "__all__"

    owner = graphene.Field(SpectrumOwnerInterface)
    color = graphene.String()

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


class FilterPlacement(DjangoObjectType):
    id = graphene.ID()
    spectrum = graphene.Field(Spectrum)
    spectrumId = id = graphene.ID()
    name = graphene.String()

    class Meta:
        model = models.FilterPlacement
        fields = "__all__"

    def resolve_spectrum(self, info):
        return self.filter.spectrum

    def resolve_name(self, info):
        return self.filter.name

    def resolve_spectrumId(self, info):
        return self.filter.spectrum.id

    def resolve_id(self, info):
        return self.filter_id


class Microscope(BaseDjangoObjectType):
    class Meta:
        model = models.Microscope
        fields = (
            "id",
            "name",
            "description",
            "extra_lights",
            "extra_cameras",
            "extra_lasers",
            "optical_configs",
        )


class OpticalConfig(DjangoObjectType):
    filters = DjangoListField(FilterPlacement, field_name="filterplacement_set")
    microscope = graphene.Field(Microscope)

    class Meta:
        model = models.OpticalConfig
        fields = (
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
