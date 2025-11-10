"""Django REST Framework mixins for ETag support."""

from __future__ import annotations

from typing import TYPE_CHECKING

from django.http import HttpResponse

from fpbase.etag_utils import generate_version_etag, parse_etag_header

if TYPE_CHECKING:
    from django.db.models import Model
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

        current_etag = generate_version_etag(*self.etag_models)
        response["ETag"] = current_etag

        if_none_match = request.headers.get("if-none-match")
        if if_none_match:
            client_etags = parse_etag_header(if_none_match)
            if "*" in client_etags or current_etag in client_etags:
                not_modified = HttpResponse(status=304)
                not_modified["ETag"] = current_etag
                for header in ("Cache-Control", "Vary", "Expires"):
                    if header in response:
                        not_modified[header] = response[header]
                return not_modified

        return response
