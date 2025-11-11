"""Django REST Framework mixins for ETag support."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fpbase.etag_utils import check_etag_match, generate_version_etag

if TYPE_CHECKING:
    from django.db.models import Model
    from django.http import HttpResponse
    from rest_framework.request import Request
    from rest_framework.response import Response


class ETagMixin:
    """Add ETag support to DRF APIView classes.

    Usage:
        class MyAPIView(ETagMixin, ListAPIView):
            etag_models = [Protein, State]
    """

    etag_models: list[type[Model]] = []

    def finalize_response(self, request: Request, response: Response, *args, **kwargs) -> Response | HttpResponse:
        response = super().finalize_response(request, response, *args, **kwargs)  # type: ignore[misc]

        if request.method not in ("GET", "HEAD") or response.status_code != 200 or not self.etag_models:
            return response

        # Add ETag header
        current_etag = generate_version_etag(*self.etag_models)
        response["ETag"] = current_etag

        # Check if client's ETag matches and return 304 if so
        not_modified = check_etag_match(request, *self.etag_models, base_response=response)
        if not_modified:
            return not_modified

        return response
