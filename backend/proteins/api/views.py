from django.db.models import F, Max, Prefetch
from django.http import HttpRequest, HttpResponse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_control, cache_page
from django.views.decorators.http import condition
from django_filters import rest_framework as filters
from rest_framework.generics import (
    ListAPIView,
    RetrieveAPIView,
    RetrieveUpdateDestroyAPIView,
)
from rest_framework.permissions import AllowAny, IsAdminUser, IsAuthenticated
from rest_framework.settings import api_settings
from rest_framework_csv import renderers as r

import proteins.models as pm
from fpbase.cache_utils import get_model_version

from ..filters import ProteinFilter, SpectrumFilter, StateFilter
from ..models.microscope import get_cached_optical_configs
from ..models.spectrum import get_cached_spectra_info
from .serializers import (
    BasicProteinSerializer,
    ProteinSerializer,
    ProteinSerializer2,
    ProteinSpectraSerializer,
    ProteinTableSerializer,
    SpectrumSerializer,
    StateSerializer,
)


def _spectra_etag(request: HttpRequest) -> str:
    """Compute weak ETag for spectra list based on model versions."""
    version = get_model_version(pm.Camera, pm.Dye, pm.Filter, pm.Light, pm.Protein, pm.Spectrum, pm.State)
    return f'W/"{version}"'


def _optical_configs_etag(request: HttpRequest) -> str:
    """Compute weak ETag for optical configs list based on model versions."""
    version = get_model_version(pm.Microscope, pm.OpticalConfig)
    return f'W/"{version}"'


@condition(etag_func=_spectra_etag)
@cache_control(public=True, max_age=300, must_revalidate=True)
def spectra_list(request: HttpRequest) -> HttpResponse:
    """Return cached spectra list with ETag support."""
    data = get_cached_spectra_info()
    return HttpResponse(
        data,
        content_type="application/json",
        headers={"Vary": "Accept-Encoding"},
    )


@condition(etag_func=_optical_configs_etag)
@cache_control(public=True, max_age=300, must_revalidate=True)
def optical_configs_list(request: HttpRequest) -> HttpResponse:
    """Return cached optical configs list with ETag support."""
    data = get_cached_optical_configs()
    return HttpResponse(
        data,
        content_type="application/json",
        headers={"Vary": "Accept-Encoding"},
    )


class SpectrumList(ListAPIView):
    queryset = pm.Spectrum.objects.all()
    serializer_class = SpectrumSerializer
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = SpectrumFilter


class SpectrumDetail(RetrieveAPIView):
    queryset = pm.Spectrum.objects.prefetch_related("owner_state")
    permission_classes = (AllowAny,)
    serializer_class = SpectrumSerializer


class ProteinListAPIView2(ListAPIView):
    queryset = pm.Protein.objects.all().prefetch_related("states", "transitions")
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer2
    lookup_field = "slug"  # Don't use Protein.id!
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = ProteinFilter
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class ProteinListAPIView(ListAPIView):
    queryset = (
        pm.Protein.objects.all()
        .prefetch_related(
            "states__spectra",  # Prefetch spectra for each state to avoid N+1 queries
            Prefetch(
                "transitions",
                queryset=pm.StateTransition.objects.select_related("from_state", "to_state"),
            ),
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
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class BasicProteinListAPIView(ProteinListAPIView):
    queryset = (
        pm.Protein.visible.filter(switch_type=pm.Protein.BASIC)
        .select_related("default_state")
        .annotate(rate=Max(F("default_state__bleach_measurements__rate")))
    )
    permission_classes = (AllowAny,)
    serializer_class = BasicProteinSerializer


class ProteinRetrieveUpdateDestroyAPIView(RetrieveUpdateDestroyAPIView):
    queryset = pm.Protein.objects.all()
    permission_classes = (IsAdminUser,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id


class ProteinRetrieveAPIView(RetrieveAPIView):
    queryset = pm.Protein.objects.all()
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id


class StatesListAPIView(ListAPIView):
    queryset = pm.State.objects.all().select_related("protein")
    permission_classes = (IsAuthenticated,)
    serializer_class = StateSerializer
    lookup_field = "slug"  # Don't use State.id!
    renderer_classes = [r.CSVRenderer, *api_settings.DEFAULT_RENDERER_CLASSES]  # pyright: ignore[reportAssignmentType]
    filter_backends = (filters.DjangoFilterBackend,)
    filterset_class = StateFilter


class ProteinSpectraListAPIView(ListAPIView):
    permission_classes = (AllowAny,)
    serializer_class = ProteinSpectraSerializer
    queryset = pm.Protein.objects.with_spectra().prefetch_related("states")


class ProteinTableAPIView(ListAPIView):
    """Optimized API endpoint for the protein table view.

    Includes efficient queries with prefetch_related and select_related to
    avoid N+1 query problems. Only returns visible proteins with their states.
    """

    queryset = (
        pm.Protein.visible.all()
        .prefetch_related(
            Prefetch("states", queryset=pm.State.objects.filter(is_dark=False))
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
