import graphene
import graphene_django_optimizer as gdo
from django.db.models import Prefetch
from graphene_django.types import DjangoObjectType

from . import models


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


class Protein(gdo.OptimizedDjangoObjectType):
    id = graphene.String()
    
    class Meta:
        model = models.Protein
        exclude_fields = ("id", "status", "status_changed", "uuid", "base_name")

    def resolve_id(self, info):
        return self.uuid


class FluorophoreInterface(graphene.Interface):
    qy = graphene.Float()
    extCoeff = graphene.Float()
    twopPeakgm = graphene.Float()
    exMax = graphene.Float()
    emMax = graphene.Float()

    def resolve_extCoeff(self, info):
        return self.ext_coeff

    def resolve_twopPeakgm(self, info):
        return self.twop_peakGM

    def resolve_exMax(self, info):
        return self.ex_max

    def resolve_emMax(self, info):
        return self.em_max


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
    class Meta:
        interfaces = (SpectrumOwnerInterface, FluorophoreInterface)
        model = models.State
        exclude_fields = ()

    protein = graphene.Field(Protein)

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
        exclude_fields = (
            "owner_state",
            "owner_dye",
            "owner_filter",
            "owner_camera",
            "owner_light",
        )

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
    path = graphene.String()
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
    state = graphene.Field(State, id=graphene.Int())
    spectrum = graphene.Field(Spectrum, id=graphene.Int())
    states = graphene.List(State)
    # spectra = graphene.List(Spectrum)
    spectra = graphene.List(SpectrumInfo)
    opticalConfigs = graphene.List(OpticalConfig)
    opticalConfig = graphene.Field(OpticalConfig, id=graphene.Int())
    proteins = graphene.List(Protein)
    protein = graphene.Field(Protein, id=graphene.String())
    microscopes = graphene.List(Microscope)

    def resolve_microscopes(self, info, **kwargs):
        return gdo.query(models.Microscope.objects.all(), info)

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

    def resolve_state(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            return models.State.objects.get(id=id)
        return None

    def resolve_spectrum(self, info, **kwargs):
        id = kwargs.get("id")
        if id is not None:
            try:
                return gdo.query(models.Spectrum.objects.filter(id=id), info).get()
            except models.Spectrum.DoesNotExist:
                return None
        return None

    def resolve_states(self, info, **kwargs):
        return models.State.objects.all()

    def resolve_spectra(self, info, **kwargs):
        requested_fields = get_requested_fields(info)
        if "owner" in requested_fields:
            return models.Spectrum.objects.sluglist()
        return models.Spectrum.objects.all().values(*requested_fields)

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
