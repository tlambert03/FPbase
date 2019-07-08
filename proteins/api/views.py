from rest_framework.generics import (
    ListAPIView,
    RetrieveUpdateDestroyAPIView,
    RetrieveAPIView,
)
from rest_framework.permissions import IsAuthenticated, IsAdminUser, AllowAny

from ..models import Protein, State, Spectrum
from .serializers import (
    ProteinSerializer,
    ProteinSerializer2,
    SpectrumSerializer,
    BasicProteinSerializer,
    StateSerializer,
    ProteinSpectraSerializer,
)

from django.db.models import F, Max
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.http import JsonResponse
from django_filters import rest_framework as filters
from ..filters import ProteinFilter, StateFilter, SpectrumFilter
from rest_framework.settings import api_settings
from rest_framework_csv import renderers as r


class SpectrumList(ListAPIView):
    queryset = Spectrum.objects.all()
    serializer_class = SpectrumSerializer
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = SpectrumFilter


class SpectrumDetail(RetrieveAPIView):
    queryset = Spectrum.objects.prefetch_related("owner_state")
    permission_classes = (AllowAny,)
    serializer_class = SpectrumSerializer


class ProteinListAPIView2(ListAPIView):
    queryset = Protein.objects.all().prefetch_related("states", "transitions")
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer2
    lookup_field = "slug"  # Don't use Protein.id!
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = ProteinFilter
    renderer_classes = [r.CSVRenderer] + api_settings.DEFAULT_RENDERER_CLASSES

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class ProteinListAPIView(ListAPIView):
    queryset = (
        Protein.objects.all()
        .prefetch_related("states", "transitions")
        .select_related("default_state")
    )
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id!
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = ProteinFilter
    renderer_classes = [r.CSVRenderer] + api_settings.DEFAULT_RENDERER_CLASSES

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class BasicProteinListAPIView(ProteinListAPIView):
    queryset = (
        Protein.visible.filter(switch_type=Protein.BASIC)
        .select_related("default_state")
        .annotate(rate=Max(F("default_state__bleach_measurements__rate")))
    )
    permission_classes = (AllowAny,)
    serializer_class = BasicProteinSerializer


class ProteinRetrieveUpdateDestroyAPIView(RetrieveUpdateDestroyAPIView):
    queryset = Protein.objects.all()
    permission_classes = (IsAdminUser,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id


class ProteinRetrieveAPIView(RetrieveAPIView):
    queryset = Protein.objects.all()
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id


class StatesListAPIView(ListAPIView):
    queryset = State.objects.all().select_related("protein")
    permission_classes = (IsAuthenticated,)
    serializer_class = StateSerializer
    lookup_field = "slug"  # Don't use State.id!
    renderer_classes = [r.CSVRenderer] + api_settings.DEFAULT_RENDERER_CLASSES
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = StateFilter


class ProteinSpectraListAPIView(ListAPIView):
    permission_classes = (AllowAny,)
    serializer_class = ProteinSpectraSerializer
    queryset = Protein.objects.with_spectra().prefetch_related("states")


def spectraslugs(request):
    return JsonResponse(Spectrum.objects.sluglist(), safe=False)
