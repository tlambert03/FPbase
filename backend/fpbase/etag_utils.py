"""ETag generation and parsing utilities."""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING, Protocol

from django.http import HttpResponseNotModified

from fpbase.cache_utils import get_model_version

if TYPE_CHECKING:
    from django.db.models import Model
    from django.http import HttpRequest, HttpResponse
    from django.http.response import HttpResponseBase


class SupportsHeaderAccess(Protocol):
    """Protocol for objects that support dict-like header access."""

    def __contains__(self, key: object, /) -> bool: ...
    def __getitem__(self, key: str, /) -> str: ...


def generate_version_etag(*model_classes: type[Model]) -> str:
    """Generate a weak ETag from model versions."""
    return f'W/"{get_model_version(*model_classes)}"'


def generate_content_etag(content: str | bytes) -> str:
    """Generate a strong ETag from response content."""
    if isinstance(content, str):
        content = content.encode("utf-8")
    return f'"{hashlib.md5(content, usedforsecurity=False).hexdigest()}"'


def parse_etag_header(header_value: str | None) -> list[str]:
    """Parse If-None-Match header into list of ETags."""
    if not header_value:
        return []
    if header_value.strip() == "*":
        return ["*"]
    return [etag.strip() for etag in header_value.split(",") if etag.strip()]


def check_etag_match(
    request: HttpRequest,
    *models: type[Model],
    base_response: HttpResponseBase | None = None,
) -> HttpResponse | None:
    """Check if client ETag matches current version and return 304 response if so.

    Args:
        request: The HTTP request object
        *models: Model classes to generate version ETag from
        extra_headers: Optional headers to copy to the 304 response

    Returns:
        HttpResponse with 304 status if ETags match, None otherwise
    """

    if request.method not in ("GET", "HEAD"):
        return None
    if not (if_none_match := request.headers.get("if-none-match")):
        return None

    current_etag = generate_version_etag(*models)
    client_etags = parse_etag_header(if_none_match)
    if "*" not in client_etags and current_etag not in client_etags:
        return None

    # ETags match - return 304 Not Modified
    response = HttpResponseNotModified(headers={"ETag": current_etag})

    # Copy over caching-related headers if provided
    if base_response:
        for header in ("Cache-Control", "Vary", "Expires"):
            if header in base_response:
                response[header] = base_response[header]

    return response
