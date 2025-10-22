from django.db.models import F, Max, Prefetch
from django.http import HttpResponse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_control, cache_page
from django_filters import rest_framework as filters
from rest_framework.generics import (
    ListAPIView,
    RetrieveAPIView,
    RetrieveUpdateDestroyAPIView,
)
from rest_framework.permissions import AllowAny, IsAdminUser, IsAuthenticated
from rest_framework.settings import api_settings
from rest_framework_csv import renderers as r

from ..filters import ProteinFilter, SpectrumFilter, StateFilter
from ..models import Protein, State
from ..models.microscope import get_cached_optical_configs
from ..models.spectrum import Spectrum, get_cached_spectra_info
from .serializers import (
    BasicProteinSerializer,
    ProteinSerializer,
    ProteinSerializer2,
    ProteinSpectraSerializer,
    ProteinTableSerializer,
    SpectrumSerializer,
    StateSerializer,
)


def spectraslugs(request):
    spectrainfo = get_cached_spectra_info()
    return HttpResponse(spectrainfo, content_type="application/json")


def ocinfo(request):
    ocinfo = get_cached_optical_configs()
    return HttpResponse(ocinfo, content_type="application/json")


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
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class ProteinListAPIView(ListAPIView):
    queryset = (
        Protein.objects.all()
        .prefetch_related(
            "states__spectra",  # Prefetch spectra for each state to avoid N+1 queries
            "transitions",
        )
        .select_related(
            "default_state",
            "primary_reference",  # Needed for DOI field in serializer
        )
    )
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id!
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = ProteinFilter
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]

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
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = StateFilter


class ProteinSpectraListAPIView(ListAPIView):
    permission_classes = (AllowAny,)
    serializer_class = ProteinSpectraSerializer
    queryset = Protein.objects.with_spectra().prefetch_related("states")


class ProteinTableAPIView(ListAPIView):
    """Optimized API endpoint for the protein table view.

    Includes efficient queries with prefetch_related and select_related to
    avoid N+1 query problems. Only returns visible proteins with their states.
    """

    queryset = (
        Protein.visible.all()
        .prefetch_related(
            Prefetch("states", queryset=State.objects.filter(is_dark=False))
        )  # Prefetch only non-dark states
        .select_related("primary_reference")  # Needed for year field
        .order_by("name")
    )
    permission_classes = (AllowAny,)
    serializer_class = ProteinTableSerializer
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = ProteinFilter

    @method_decorator(cache_control(public=True, max_age=600))
    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)
