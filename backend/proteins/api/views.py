import difflib

from django.db.models import F, Max, Prefetch, Q
from django.http import Http404, HttpRequest, HttpResponse, JsonResponse
from django.urls import reverse
from django.utils.cache import get_conditional_response
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_control
from django.views.decorators.http import condition
from django_filters import rest_framework as filters
from drf_spectacular.utils import OpenApiParameter, extend_schema, extend_schema_view
from rest_framework.exceptions import NotFound, ValidationError
from rest_framework.generics import (
    ListAPIView,
    RetrieveAPIView,
    RetrieveUpdateDestroyAPIView,
)
from rest_framework.pagination import LimitOffsetPagination, _positive_int
from rest_framework.permissions import AllowAny, IsAdminUser, IsAuthenticated
from rest_framework.response import Response
from rest_framework.settings import api_settings
from rest_framework_csv import renderers as r

import proteins.models as pm
from fpbase.cache_utils import cache_page_by_data_version, get_model_version
from fpbase.views import ExpensiveListAnonThrottle, SameOriginExemptAnonThrottle, UserThrottle
from proteins.api._tweaks import ModelSerializer
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
# the API docs tell users to copy; `fields`/`include_fields` choose the fields returned.
NON_FILTER_PARAMS = {api_settings.URL_FORMAT_OVERRIDE, "display", "fields", "include_fields"}


class StrictDjangoFilterBackend(filters.DjangoFilterBackend):
    """Reject unknown query params, instead of ignoring them.

    An ignored param (`?search=`, `?page=`) means an unfiltered response: the client
    guessing at our API gets the whole database, every time.
    """

    def filter_queryset(self, request, queryset, view):
        allowed = {*self.get_filterset_class(view, queryset).base_filters, *NON_FILTER_PARAMS}
        if paginator := view.paginator:
            allowed |= set(paginator.query_params)
        if not issubclass(view.get_serializer_class(), ModelSerializer):
            # a plain serializer would ignore them: say so, rather than answer in full
            allowed -= {"fields", "include_fields"}
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

    @method_decorator(cache_page_by_data_version())
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class OptionalLimitOffsetPagination(LimitOffsetPagination):
    """Honor ?limit=&offset= (or ?page=&page_size=) when given, keeping a bare list.

    Without either the full list is returned, as it always has been.  Clients that
    page until they receive an empty list will now actually terminate.
    """

    default_limit = None
    max_limit = 1000
    # the other common spelling: page=2&page_size=50 is limit=50&offset=50
    page_query_param = "page"
    page_size_query_param = "page_size"
    default_page_size = 100

    @property
    def query_params(self) -> tuple[str, ...]:
        return (
            self.limit_query_param,
            self.offset_query_param,
            self.page_query_param,
            self.page_size_query_param,
        )

    def get_limit(self, request):
        params = request.query_params
        if self.limit_query_param not in params and (
            self.page_size_query_param in params or self.page_query_param in params
        ):
            try:
                return _positive_int(
                    params.get(self.page_size_query_param, self.default_page_size),
                    strict=True,
                    cutoff=self.max_limit,
                )
            except (KeyError, ValueError):
                return self.default_page_size
        return super().get_limit(request)

    def get_offset(self, request):
        params = request.query_params
        if self.offset_query_param not in params and self.page_query_param in params:
            try:
                page = _positive_int(params[self.page_query_param], strict=True)
            except (KeyError, ValueError):
                page = 1
            return (page - 1) * (self.limit or 0)
        return super().get_offset(request)

    def get_paginated_response(self, data):
        return Response(data)

    def get_schema_operation_parameters(self, view):
        return [
            *super().get_schema_operation_parameters(view),
            {
                "name": self.page_query_param,
                "required": False,
                "in": "query",
                "description": "Page number (1-based); the same as offset=(page-1)*page_size.",
                "schema": {"type": "integer"},
            },
            {
                "name": self.page_size_query_param,
                "required": False,
                "in": "query",
                "description": f"Results per page (default {self.default_page_size}).",
                "schema": {"type": "integer"},
            },
        ]


FIELDS_PARAMETERS = [
    OpenApiParameter(
        "fields",
        str,
        description="Only these fields, comma-separated; nested fields as `states__ex_max`.",
    ),
    OpenApiParameter(
        "include_fields",
        str,
        description="Also these on-demand fields: `states__ex_spectrum`, `states__em_spectrum`.",
    ),
]


@extend_schema_view(
    get=extend_schema(
        description=(
            "Proteins matching the query, as a list. Without `limit` or `page`, every "
            "protein is returned. See https://www.fpbase.org/api/ for the filters."
        ),
        parameters=FIELDS_PARAMETERS,
    )
)
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

    @method_decorator(cache_page_by_data_version())
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


@extend_schema_view(
    get=extend_schema(
        description="One protein.",
        parameters=[
            OpenApiParameter(
                "slug",
                str,
                OpenApiParameter.PATH,
                description=(
                    "The protein's slug (the last part of its page URL). An FPbase ID, "
                    "a name or alias, or a PDB ID is accepted too."
                ),
            ),
            *FIELDS_PARAMETERS,
        ],
    )
)
class ProteinRetrieveAPIView(RetrieveAPIView):
    queryset = ProteinListAPIView.queryset
    permission_classes = (AllowAny,)
    serializer_class = ProteinSerializer
    lookup_field = "slug"  # Don't use Protein.id

    @method_decorator(cache_page_by_data_version())
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_object(self):
        # slugs are lowercase, but clients tend to ask for `/api/proteins/mCherry/`
        key = self.kwargs["slug"]
        self.kwargs["slug"] = key.lower()
        try:
            return super().get_object()
        except Http404:
            return self._get_by_other_identifier(key)

    def _get_by_other_identifier(self, key: str) -> pm.Protein:
        """The protein a client meant by an FPbase ID, a name or alias, or a PDB ID.

        (Clients try all of these here: `/api/proteins/2IB5/`, `/api/proteins/R9NL8/`.)
        """
        if "\x00" in key:  # (postgres rejects it)
            raise NotFound(f"No protein matches {key!r}")
        queryset = self.filter_queryset(self.get_queryset())
        lookups = (
            Q(uuid__iexact=key),
            Q(name__iexact=key) | Q(aliases__icontains=key),
            Q(pdb__contains=[key.upper()]),
        )
        for lookup in lookups:
            # (a name lookup's `icontains` only narrows: the alias must match exactly)
            candidates = pm.Protein.objects.filter(lookup).only(*pm.PROTEIN_NAME_FIELDS)
            ids = [p.id for p in candidates if pm.protein_is_named(p, key)]
            matches = list(queryset.filter(id__in=ids)[:3])
            if len(matches) == 1:
                return matches[0]
            if matches:
                slugs = ", ".join(sorted(p.slug for p in matches))
                raise NotFound(
                    f"{key!r} matches more than one protein ({slugs}): "
                    f"request one by slug, or list them with ?pdb={key}"
                )
        raise NotFound(f"No protein matches {key!r} (as a slug, FPbase ID, name, alias or PDB ID)")


class StatesListAPIView(ListAPIView):
    queryset = pm.State.objects.all().select_related("protein")
    permission_classes = (IsAuthenticated,)
    serializer_class = StateSerializer
    lookup_field = "slug"  # Don't use State.id!
    renderer_classes = [*api_settings.DEFAULT_RENDERER_CLASSES, r.CSVRenderer]  # pyright: ignore[reportAssignmentType]
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = StateFilter


class SpectraPagination(OptionalLimitOffsetPagination):
    max_limit = 200


@extend_schema_view(
    get=extend_schema(
        description=(
            "Proteins' spectra. Every protein's spectra at once is more than the server "
            "can build: filter (e.g. `name=`), or page with `limit` (at most 200)."
        )
    )
)
class ProteinSpectraListAPIView(ListAPIView):
    permission_classes = (AllowAny,)
    serializer_class = ProteinSpectraSerializer
    queryset = pm.Protein.objects.with_spectra().prefetch_related("states__spectra")
    # without these, every filtered query returned (and serialized) every spectrum
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter
    pagination_class = SpectraPagination
    throttle_classes = [ExpensiveListAnonThrottle, *api_settings.DEFAULT_THROTTLE_CLASSES]  # pyright: ignore[reportAssignmentType]

    @method_decorator(cache_page_by_data_version())
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def list(self, request, *args, **kwargs):
        # the whole dump (~6 MB of JSON, far more in memory) has crashed the dyno
        paging = {*self.paginator.query_params}  # pyright: ignore[reportOptionalMemberAccess]
        filters = set(request.query_params) - NON_FILTER_PARAMS - paging
        if not filters and not (paging & set(request.query_params)):
            raise ValidationError(
                {
                    "detail": "This would be every protein's spectra at once: filter the "
                    "proteins (e.g. ?name=mCherry), or page through them with ?limit= "
                    f"(at most {SpectraPagination.max_limit}) and ?offset=.",
                    "docs": request.build_absolute_uri(reverse("api:api")),
                }
            )
        return super().list(request, *args, **kwargs)


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
    throttle_classes = [SameOriginExemptAnonThrottle, UserThrottle]
    filter_backends = (StrictDjangoFilterBackend,)
    filterset_class = ProteinAPIFilter

    @method_decorator(cache_page_by_data_version(public=True))
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)
