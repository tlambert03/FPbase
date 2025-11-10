"""ETag generation and parsing utilities."""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

from fpbase.cache_utils import get_model_version

if TYPE_CHECKING:
    from django.db.models import Model


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
