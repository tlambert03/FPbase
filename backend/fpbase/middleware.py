import logging

from django.conf import settings
from django.core.exceptions import MiddlewareNotUsed
from django.http import HttpResponseForbidden
from django.shortcuts import redirect

logger = logging.getLogger(__name__)


# https://github.com/etianen/django-herokuapp/blob/master/herokuapp/middleware.py
class CanonicalDomainMiddleware:
    """Middleware that redirects to a canonical domain."""

    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.
        if settings.DEBUG or not settings.CANONICAL_URL:
            logger.debug("Not using CANONICAL_URL")
            raise MiddlewareNotUsed

    def __call__(self, request):
        """If the request domain is not the canonical domain, redirect."""
        hostname = request.get_host().split(":", 1)[0]
        # Don't perform redirection for testing or local development.
        if hostname in ("testserver", "localhost", "127.0.0.1"):
            return self.get_response(request)
        if "googleaecf5301782589e7" in request.get_full_path():
            return self.get_response(request)
        # Check against the site domain.
        canonical_hostname = settings.CANONICAL_URL.split("://", 1)[1]
        if hostname != canonical_hostname:
            canonical_url = settings.CANONICAL_URL + request.get_full_path()
            logger.debug(f"REDIRECTING {request.get_host()} to {canonical_url}")
            return redirect(canonical_url, permanent=True)
        return self.get_response(request)


class BlackListMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response
        self.blacklist: set[str] = settings.BLOCKED_IPS

    def __call__(self, request):
        if request.META.get("REMOTE_ADDR") in self.blacklist:
            return HttpResponseForbidden()

        return self.get_response(request)
