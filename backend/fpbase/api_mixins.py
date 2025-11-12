"""Django REST Framework mixins for ETag support."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar

from rest_framework.views import APIView

from fpbase.etag_utils import check_etag_match, etagged_response

if TYPE_CHECKING:
    from collections.abc import Sequence

    from django.db.models import Model
    from django.http import HttpResponse
    from rest_framework.request import Request
    from rest_framework.response import Response as DRFResponse


class ETagMixin(APIView):
    """Add ETag support to DRF APIView classes.

    Usage:
        class MyAPIView(ETagMixin, ListAPIView):
            etag_models = [Protein, State]
    """

    etag_models: ClassVar[Sequence[type[Model]]] = ()

    def finalize_response(
        self, request: Request, response: DRFResponse, *args: Any, **kwargs: Any
    ) -> DRFResponse | HttpResponse:
        response = super().finalize_response(request, response, *args, **kwargs)
        # Check if client's ETag matches and return 304 if so
        if (
            response.status_code == 200
            and request.method in ("GET", "HEAD")
            and self.etag_models
            and (not_modified := check_etag_match(request, *self.etag_models, base_response=response))
        ):
            return not_modified

        # Otherwise, add ETag to 200 response
        return etagged_response(response, request, *self.etag_models)
