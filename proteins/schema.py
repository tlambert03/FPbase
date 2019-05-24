import graphene
import graphene_django_optimizer as gdo
from graphene_django.types import DjangoObjectType


from . import models


class SpectrumOwnerInterface(graphene.Interface):
    name = graphene.String()
    id = graphene.ID()
    typ = graphene.String()
    slug = graphene.String()

    def resolve_typ(self, info):
        return self.__class__.__name__.lower()


class Camera(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
        model = models.Camera


class Dye(DjangoObjectType):
    class Meta:
        interfaces = (SpectrumOwnerInterface,)
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
        interfaces = (SpectrumOwnerInterface,)
        model = models.State

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


class Query(graphene.ObjectType):
    state = graphene.Field(State, id=graphene.Int())
    spectrum = graphene.Field(Spectrum, id=graphene.Int())
    states = graphene.List(State)
    spectra = graphene.List(Spectrum)

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
        return gdo.query(models.Spectrum.objects.all(), info)

