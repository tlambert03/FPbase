import logging

from django.conf import settings
from django.contrib.admin.views.decorators import staff_member_required
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.views.generic.edit import FormView
from graphene_django.views import GraphQLView
from rest_framework import exceptions
from rest_framework.settings import api_settings
from sentry_sdk import last_event_id

from fpbase.forms import ContactForm
from proteins.models import Protein, Spectrum

logger = logging.getLogger(__name__)


class RateLimitedGraphQLView(GraphQLView):
    """
    GraphQL view with rate limiting using DRF's throttle infrastructure.

    Leverages Django REST Framework's battle-tested throttling system:
    - Uses DEFAULT_THROTTLE_CLASSES from settings (AnonRateThrottle, UserRateThrottle)
    - Automatically handles X-Forwarded-For for Heroku deployments
    - Raises DRF's Throttled exception which includes retry-after information
    - Converts the exception to GraphQL error format with proper HTTP headers
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
            retry_after = int(exc.wait) if exc.wait else 60  # ty: ignore[unresolved-attribute]

            # Log rate limit event with structured data
            logger.warning(
                "GraphQL rate limit exceeded",
                extra={
                    "user_id": request.user.id if request.user.is_authenticated else None,
                    "is_authenticated": request.user.is_authenticated,
                    "path": request.path,
                    "method": request.method,
                    "user_agent": request.headers.get("user-agent", "")[:200],
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


@cache_page(10)
def testview(request):
    import logging

    logger = logging.getLogger(__name__)
    p = Protein.objects.get(name="mNeonGreen")
    logger.info(p)
    return render(request, "pages/test.html", {"protein": p})
