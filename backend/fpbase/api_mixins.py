"""Django REST Framework mixins for ETag support.

This module provides mixins that add HTTP conditional request support
(ETags, If-None-Match, 304 Not Modified) to DRF API views.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from django.http import HttpResponse

from fpbase.etag_utils import generate_version_etag, parse_etag_header

if TYPE_CHECKING:
    from django.db.models import Model
    from rest_framework.request import Request
    from rest_framework.response import Response


class ETagMixin:
    """Add ETag support to Django REST Framework APIView classes.

    This mixin adds automatic ETag generation based on model versions and
    handles If-None-Match conditional requests, returning 304 Not Modified
    when appropriate.

    Usage
    -----
    Simply add this mixin to your APIView and specify which models to track:

    >>> class ProteinListAPIView(ETagMixin, ListAPIView):
    ...     etag_models = [Protein, State]
    ...     # ... rest of your view code ...

    When to use 304 responses
    -------------------------
    - GET requests with matching If-None-Match header
    - HEAD requests with matching If-None-Match header
    - Not for POST/PUT/DELETE (those modify resources)

    Attributes
    ----------
    etag_models : list[type[Model]]
        List of Django model classes to track for version-based ETags.
        The ETag will be generated from the combined versions of all models.
        Default: [] (no models tracked, no ETag generated)
    """

    etag_models: list[type[Model]] = []

    def finalize_response(self, request: Request, response: Response, *args, **kwargs) -> Response | HttpResponse:
        """Add ETag header and handle conditional requests.

        This is called by DRF after the view has generated the response
        but before it's sent to the client.

        Parameters
        ----------
        request
            The DRF request object.
        response
            The DRF response object.
        *args, **kwargs
            Additional arguments passed by DRF.

        Returns
        -------
        Response | HttpResponse
            Either the original response (200) or a 304 Not Modified response.
        """
        # Call parent first to ensure response is finalized
        response = super().finalize_response(request, response, *args, **kwargs)  # type: ignore[misc]

        # Only add ETags for successful GET/HEAD requests
        if request.method not in ("GET", "HEAD") or response.status_code != 200:
            return response

        # Only add ETags if models are specified
        if not self.etag_models:
            return response

        # Generate ETag from model versions
        current_etag = generate_version_etag(*self.etag_models)

        # Add ETag header to response
        response["ETag"] = current_etag

        # Check if client sent If-None-Match header
        if_none_match = request.headers.get("if-none-match")
        if if_none_match:
            client_etags = parse_etag_header(if_none_match)

            # Check if current ETag matches any client ETags
            # or if client sent wildcard (*)
            if "*" in client_etags or current_etag in client_etags:
                # Data hasn't changed, return 304 Not Modified
                # Per RFC 7232, 304 response should include the ETag header
                not_modified = HttpResponse(status=304)
                not_modified["ETag"] = current_etag

                # Copy cache-related headers from original response
                for header in ("Cache-Control", "Vary", "Expires"):
                    if header in response:
                        not_modified[header] = response[header]

                return not_modified

        return response
