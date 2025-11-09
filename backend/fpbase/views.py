import logging

from django.conf import settings
from django.contrib.admin.views.decorators import staff_member_required
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django.views.generic.edit import FormView
from graphene_django.views import GraphQLView
from rest_framework import exceptions
from rest_framework.settings import api_settings
from rest_framework.throttling import AnonRateThrottle
from sentry_sdk import last_event_id

from fpbase.etag_utils import generate_version_etag, parse_etag_header
from fpbase.forms import ContactForm
from proteins.models import OpticalConfig, Protein, Spectrum

logger = logging.getLogger(__name__)


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

    Features:
    - Rate limiting using DRF's throttle infrastructure
    - ETag-based caching for conditional requests (304 Not Modified)
    - Automatic cache invalidation when Spectrum or OpticalConfig changes

    The ETag is based on Spectrum and OpticalConfig model versions, since
    these are the primary data returned by the spectra viewer GraphQL queries.
    """

    # Use the same throttle classes as the REST API (from settings.REST_FRAMEWORK)
    throttle_classes = api_settings.DEFAULT_THROTTLE_CLASSES

    # Models to track for ETag generation
    etag_models = [Spectrum, OpticalConfig]

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

        # Check ETag for conditional requests (304 Not Modified)
        # Only for GET/POST requests (GraphQL uses POST for queries)
        if request.method in ("GET", "POST"):
            current_etag = generate_version_etag(*self.etag_models)
            if_none_match = request.headers.get("if-none-match")

            if if_none_match:
                client_etags = parse_etag_header(if_none_match)
                if "*" in client_etags or current_etag in client_etags:
                    # Data hasn't changed, return 304 Not Modified
                    response = HttpResponse(status=304)
                    response["ETag"] = current_etag
                    return response

        # Process the GraphQL request
        response = super().dispatch(request, *args, **kwargs)

        # Add ETag header to successful responses
        if response.status_code == 200 and request.method in ("GET", "POST"):
            current_etag = generate_version_etag(*self.etag_models)
            response["ETag"] = current_etag

        return response


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
