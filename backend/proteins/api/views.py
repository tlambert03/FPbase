import difflib

from django.db.models import F, Max, Prefetch
from django.http import HttpRequest, HttpResponse, JsonResponse
from django.urls import reverse
from django.utils.cache import get_conditional_response
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_control, cache_page
from django.views.decorators.http import condition
from django_filters import rest_framework as filters
from rest_framework.exceptions import ValidationError
from rest_framework.generics import (
    ListAPIView,
    RetrieveAPIView,
    RetrieveUpdateDestroyAPIView,
)
from rest_framework.pagination import LimitOffsetPagination
from rest_framework.permissions import AllowAny, IsAdminUser, IsAuthenticated
from rest_framework.response import Response
from rest_framework.settings import api_settings
from rest_framework.throttling import UserRateThrottle
from rest_framework_csv import renderers as r

import proteins.models as pm
from fpbase.cache_utils import get_model_version
from fpbase.views import ExpensiveListAnonThrottle, SameOriginExemptAnonThrottle
from proteins.api.serializers import (
    BasicProteinSerializer,
    ProteinSerializer,
    ProteinSerializer2,
    ProteinSpectraSerializer,
    ProteinTableSerializer,
    SpectrumSerializer,
    StateSerializer,
)
from proteins.filters import ProteinAPIFilter, SpectrumFilter, StateFilter
from proteins.models.microscope import get_cached_optical_configs
from proteins.models.spectrum import get_cached_spectra_info
from proteins.search_index import get_search_index


def _spectra_etag(request: HttpRequest) -> str:
    """Compute weak ETag for spectra list based on model versions."""
    version = get_model_version(
        pm.Camera, pm.Dye, pm.Filter, pm.Light, pm.Protein, pm.Spectrum, pm.State
    )
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


@cache_control(public=True, max_age=60 * 60)
def search_index(request: HttpRequest) -> HttpResponse:
    """Return the client-side autocomplete search index with ETag support."""
    data, etag = get_search_index()
    if response := get_conditional_response(request, etag=etag):
        return response
    return HttpResponse(data, content_type="application/json", headers={"ETag": etag})


def api_not_found(request: HttpRequest) -> JsonResponse:
    """JSON 404 for unknown API paths, pointing scripts at the real API."""
    proteins = request.build_absolute_uri(reverse("api:protein-api"))
    return JsonResponse(
        {
            "detail": f"No API endpoint at {request.path}",
            "docs": request.build_absolute_uri(reverse("api:api")),
            "proteins": proteins,
            "examples": [
                f"{proteins}mcherry/?format=json",
                f"{proteins}?name=mCherry&format=json",
            ],
            "graphql": request.build_absolute_uri("/graphql/"),
        },
        status=404,
    )


# query params that are not filters.  `display` is added to search page URLs, which
# the API docs tell users to copy.
NON_FILTER_PARAMS = {api_settings.URL_FORMAT_OVERRIDE, "display"}


class StrictDjangoFilterBackend(filters.DjangoFilterBackend):
    """Reject unknown query params, instead of ignoring them.

    An ignored param (`?search=`, `?page=`) means an unfiltered response: the client
    guessing at our API gets the whole database, every time.
    """

    def filter_queryset(self, request, queryset, view):
        allowed = {*self.get_filterset_class(view, queryset).base_filters, *NON_FILTER_PARAMS}
        if paginator := view.paginator:
            allowed |= {paginator.limit_query_param, paginator.offset_query_param}
        if unknown := sorted(set(request.query_params) - allowed):
            detail = f"Unknown query parameter(s): {', '.join(unknown)}"
            hints = {p: s for p in unknown if (s := _suggest_param(p, allowed - {"display"}))}
            if hints:
                detail += ". Did you mean: " + ", ".join(
                    f"{s} (for {p})" for p, s in hints.items()
                )
            raise ValidationError(
                {
                    "detail": detail,
                    "did_you_mean": hints,
                    "valid_parameters": sorted(allowed - {"display"}),
                    "docs": request.build_absolute_uri(reverse("api:api")),
                }
            )
        return super().filter_queryset(request, queryset, view)


def _suggest_param(param: str, allowed: set[str]) -> str | None:
    """The valid param a guess most likely meant, e.g. `pdb_id` -> `pdb`."""
    # a param on a related model: `ex_max__gte` -> `default_state__ex_max__gte`
    if suffixed := sorted(a for a in allowed if a.endswith(f"__{param}")):
        return suffixed[0]
    # an explicit exact lookup, where the bare name is exact: `slug__iexact` -> `slug`
    field, _, lookup = param.rpartition("__")
    if lookup in ("exact", "iexact") and field in allowed:
        return field
    # same field, other spelling: `pdb_id` -> `pdb`, `name__contains` -> `name__icontains`
    root = param.split("__")[0].removesuffix("_id")
    related = [a for a in allowed if a == root or a.startswith(f"{root}__")]
    # otherwise only a near-typo (`ex_maxx`), so an unrelated guess gets no hint
    matches = difflib.get_close_matches(param, related, n=1, cutoff=0)
    return (matches or difflib.get_close_matches(param, allowed, n=1, cutoff=0.8) or [None])[0]


class SpectrumList(ListAPIView):
    queryset = pm.Spectrum.objects.all()
    serializer_class = SpectrumSerializer
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = SpectrumFilter


class SpectrumDetail(RetrieveAPIView):
    queryset = pm.Spectrum.objects.prefetch_related("owner_fluor")
    permission_classes = (AllowAny,)
    serializer_class = SpectrumSerializer


class ProteinListAPIView2(ListAPIView):
    queryset = pm.Protein.objects.all().prefetch_related("states", "transitions")
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer2
    lookup_field = "slug"  # Don't use Protein.id!
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter
    renderer_classes = [*api_settings.DEFAULT_RENDERER_CLASSES, r.CSVRenderer]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class OptionalLimitOffsetPagination(LimitOffsetPagination):
    """Honor ?limit=&offset= when given, keeping the response a bare list.

    Without `limit` the full list is returned, as it always has been.  Clients that
    page until they receive an empty list will now actually terminate.
    """

    default_limit = None
    max_limit = 1000

    def get_paginated_response(self, data):
        return Response(data)


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
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter
    pagination_class = OptionalLimitOffsetPagination
    throttle_classes = [ExpensiveListAnonThrottle, *api_settings.DEFAULT_THROTTLE_CLASSES]  # pyright: ignore[reportAssignmentType]
    renderer_classes = [*api_settings.DEFAULT_RENDERER_CLASSES, r.CSVRenderer]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class BasicProteinListAPIView(ProteinListAPIView):
    queryset = (
        pm.Protein.visible.filter(switch_type=pm.Protein.SwitchingChoices.BASIC)
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
    queryset = ProteinListAPIView.queryset
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_object(self):
        # slugs are lowercase, but clients tend to ask for `/api/proteins/mCherry/`
        self.kwargs["slug"] = self.kwargs["slug"].lower()
        return super().get_object()


class StatesListAPIView(ListAPIView):
    queryset = pm.State.objects.all().select_related("protein")
    permission_classes = (IsAuthenticated,)
    serializer_class = StateSerializer
    lookup_field = "slug"  # Don't use State.id!
    renderer_classes = [*api_settings.DEFAULT_RENDERER_CLASSES, r.CSVRenderer]  # pyright: ignore[reportAssignmentType]
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = StateFilter


class ProteinSpectraListAPIView(ListAPIView):
    permission_classes = (AllowAny,)
    serializer_class = ProteinSpectraSerializer
    queryset = pm.Protein.objects.with_spectra().prefetch_related("states__spectra")
    # without these, every filtered query returned (and serialized) every spectrum
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter
    pagination_class = OptionalLimitOffsetPagination
    throttle_classes = [ExpensiveListAnonThrottle, *api_settings.DEFAULT_THROTTLE_CLASSES]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


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
    # fetched by the table page itself
    throttle_classes = [SameOriginExemptAnonThrottle, UserRateThrottle]
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter

    @method_decorator(cache_control(public=True, max_age=600))
    @method_decorator(cache_page(60 * 10))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)
