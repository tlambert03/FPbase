import hashlib
import logging

from django.conf import settings
from django.contrib.admin.views.decorators import staff_member_required
from django.http import HttpResponseNotModified, JsonResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django.views.generic.edit import FormView
from graphene_django.views import GraphQLView
from rest_framework import exceptions
from rest_framework.settings import api_settings
from rest_framework.throttling import AnonRateThrottle
from sentry_sdk import last_event_id

from fpbase.cache_utils import get_model_version
from fpbase.etag_utils import parse_etag_header
from fpbase.forms import ContactForm
from proteins.models import (
    Camera,
    Dye,
    Filter,
    Light,
    Microscope,
    OpticalConfig,
    Protein,
    Spectrum,
    State,
)

logger = logging.getLogger(__name__)


# GraphQL operation names to Django models mapping for ETag generation
# IMPORTANT: These operation names must match the query names in packages/spectra/src/api/queries.ts
# See test_graphql_operation_names_match_frontend for validation
GRAPHQL_OPERATION_ETAG_MODELS = {
    # SpectraList query: Returns list of all spectra with owner info
    # Owner can be: State (protein), Dye, Camera, Light, or Filter
    "_FPB_SpectraList": [Spectrum, State, Protein, Dye, Camera, Light, Filter],
    # OpticalConfigList query: Returns list of optical configs with microscope info
    "_FPB_OpticalConfigList": [OpticalConfig, Microscope],
    # Single spectrum query with owner info
    "_FPB_Spectrum": [Spectrum, State, Protein, Dye, Camera, Light, Filter],
    # Batch spectra queries
    "_FPB_BatchSpectra": [Spectrum, State, Protein, Dye, Camera, Light, Filter],
    # Single optical config with filters, light, camera, microscope
    "_FPB_OpticalConfig": [OpticalConfig, Microscope, Filter, Light, Camera, Spectrum],
}


def normalize_graphql_query(query: str) -> str:
    """Normalize GraphQL query for consistent hashing.

    Removes whitespace variations so semantically identical queries hash the same.
    """
    return query.strip().replace("\n", " ").replace("\t", " ").replace("  ", " ")


def generate_graphql_etag(query: str, *model_classes: type) -> str:
    """Generate ETag from GraphQL query content + model versions.

    Format: "query_hash-model_version"
    - query_hash: First 8 chars of MD5 hash of normalized query
    - model_version: Combined version hash of all models

    This ensures:
    - Different queries get different ETags (safety)
    - Same query with unchanged models reuses cache (efficiency)
    """
    normalized = normalize_graphql_query(query)
    query_hash = hashlib.md5(normalized.encode(), usedforsecurity=False).hexdigest()[:8]
    model_version = get_model_version(*model_classes)
    return f'W/"{query_hash}-{model_version}"'


class SameOriginExemptAnonThrottle(AnonRateThrottle):
    """
    Throttle class that exempts same-origin requests from rate limiting.

    This allows the FPbase spectra viewer (https://www.fpbase.org/spectra/)
    to make unlimited GraphQL requests to its own backend, while still
    throttling external API consumers.

    Same-origin is determined by checking if the Referer header matches
    the request host.
    """

    def allow_request(self, request, view):
        """Check if request should be throttled, exempting same-origin requests."""
        # Check if this is a same-origin request by comparing the referer with the host
        referer = request.headers.get("referer", "")

        # If the referer contains our host (accounting for port differences), it's a same-origin request
        # Note: This checks for the host in the referer URL (e.g., "https://www.fpbase.org/...")
        # We strip the port from host comparison to handle localhost:8000 vs fpbase.org
        if referer:
            host = request.get_host()
            # Extract the main host without port for comparison
            host_without_port = host.split(":")[0]
            # Check if host (with or without port) appears in the referer
            if host in referer or host_without_port in referer:
                return True

        # For all other requests, apply normal throttling
        return super().allow_request(request, view)


class RateLimitedGraphQLView(GraphQLView):
    """GraphQL view with rate limiting and ETag support.

    Leverages Django REST Framework's battle-tested throttling system:
    - Uses DEFAULT_THROTTLE_CLASSES from settings (AnonRateThrottle, UserRateThrottle)
    - Automatically handles X-Forwarded-For for Heroku deployments
    - Raises DRF's Throttled exception which includes retry-after information
    - Converts the exception to GraphQL error format with proper HTTP headers

    ETag support for specific GraphQL queries:
    - Uses operation name + query hash + model versions for accurate caching
    - Only applies to whitelisted operations (see GRAPHQL_OPERATION_ETAG_MODELS)
    - Safe: Different query bodies always get different ETags
    - Efficient: Only regenerates when query changes OR models change
    """

    # Use the same throttle classes as the REST API (from settings.REST_FRAMEWORK)
    throttle_classes = api_settings.DEFAULT_THROTTLE_CLASSES

    def get_throttles(self):
        """Instantiate and return the list of throttles that this view uses."""
        throttles = []
        for throttle_class in self.throttle_classes:
            try:
                throttles.append(throttle_class())
            except Exception as e:
                # If a throttle class is improperly configured, skip it
                logger.error("Error instantiating throttle %s: %s", throttle_class, str(e))
        return throttles

    def check_throttles(self, request):
        """
        Check if request should be throttled.
        Raises exceptions.Throttled if the request is throttled.

        This is adapted from rest_framework.views.APIView.check_throttles()
        """
        throttle_durations = []
        for throttle in self.get_throttles():
            if not throttle.allow_request(request, self):
                throttle_durations.append(throttle.wait())

        if throttle_durations:
            # Filter out None values (can happen with config changes)
            durations = [duration for duration in throttle_durations if duration is not None]
            duration = max(durations, default=None)

            # Raise DRF's Throttled exception (includes wait time)
            raise exceptions.Throttled(wait=duration)

    def dispatch(self, request, *args, **kwargs):
        try:
            # Check rate limits using DRF's infrastructure
            self.check_throttles(request)
        except exceptions.Throttled as exc:
            # Extract wait time from DRF's exception
            retry_after = int(exc.wait) if exc.wait else 60

            # Log rate limit event with structured data
            logger.warning(
                "GraphQL rate limit exceeded",
                extra={
                    "user_id": request.user.id if request.user.is_authenticated else None,
                    "is_authenticated": request.user.is_authenticated,
                    "path": request.path,
                    "method": request.method,
                    "user_agent": request.headers.get("user-agent", "")[:200],
                    "referer": request.headers.get("referer", "")[:200],  # for analysis
                    "retry_after": retry_after,
                    "exception_detail": str(exc.detail),
                },
            )

            # Create GraphQL-formatted error response
            response = JsonResponse(
                {
                    "errors": [
                        {
                            "message": str(exc.detail),
                            "extensions": {
                                "code": "RATE_LIMIT_EXCEEDED",
                                "retryAfter": retry_after,
                            },
                        }
                    ]
                },
                status=exc.status_code,  # 429 from DRF's Throttled exception
            )

            # Add Retry-After header (DRF's exception handler would do this too)
            response["Retry-After"] = str(retry_after)

            return response

        # ETag support for specific GraphQL GET queries
        if request.method == "GET":
            operation_name = request.GET.get("operationName")
            query = request.GET.get("query")

            # Check if this is a whitelisted operation with ETag support
            if operation_name and query and (etag_models := GRAPHQL_OPERATION_ETAG_MODELS.get(operation_name)):
                # Generate ETag from query hash + model versions
                current_etag = generate_graphql_etag(query, *etag_models)

                # Check if client's ETag matches - return 304 if so
                if if_none_match := request.headers.get("if-none-match"):
                    client_etags = parse_etag_header(if_none_match)
                    if current_etag in client_etags:
                        response = HttpResponseNotModified()
                        response["ETag"] = current_etag
                        response["Cache-Control"] = "public, max-age=600"
                        response["Vary"] = "Accept-Encoding, Origin"
                        return response

                # Execute query and add ETag to response
                response = super().dispatch(request, *args, **kwargs)
                if response.status_code == 200:
                    response["ETag"] = current_etag
                    response["Cache-Control"] = "public, max-age=600"
                    response["special-note"] = "added for debugging"
                    response["Vary"] = "Accept-Encoding, Origin"
                return response

        # No ETag support for this request (POST, unknown operation, etc.)
        return super().dispatch(request, *args, **kwargs)


class HomeView(TemplateView):
    template_name = "pages/home.html"

    def get_context_data(self):
        data = super().get_context_data()
        data["stats"] = {
            "proteins": Protein.objects.count(),
            "protspectra": Spectrum.objects.exclude(owner_state=None).count(),
        }
        return data


class ContactView(FormView):
    template_name = "pages/contact.html"
    form_class = ContactForm
    success_url = "/thanks/"

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        form.send_email()
        return super().form_valid(form)


@staff_member_required
def test500(request):
    # Return an "Internal Server Error" 500 response code.
    raise Exception("Make response code 500!")


def server_error(request, *args, **argv):
    return render(
        request,
        "500.html",
        {
            "sentry_event_id": last_event_id(),
            "sentry_dsn": getattr(settings, "SENTRY_DSN", ""),
        },
        status=500,
    )
