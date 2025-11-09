"""ETag generation and parsing utilities.

This module provides utilities for generating and parsing HTTP ETags for
conditional requests (If-None-Match, 304 Not Modified).
"""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

from fpbase.cache_utils import get_model_version

if TYPE_CHECKING:
    from django.db.models import Model


def generate_version_etag(*model_classes: type[Model]) -> str:
    """Generate a weak ETag from model versions.

    Weak ETags (W/"...") indicate that the resources are semantically equivalent
    but not byte-for-byte identical. This is perfect for version-based caching
    where the data hasn't changed, even if the JSON representation might differ
    slightly (whitespace, field order, etc.).

    Parameters
    ----------
    *model_classes
        One or more Django model classes to generate ETag from.

    Returns
    -------
    str
        Weak ETag in format: W/"<hash>"

    Examples
    --------
    >>> generate_version_etag(Spectrum)
    'W/"a3c2f1e8b9d4..."'
    >>> generate_version_etag(Spectrum, Protein)
    'W/"b4d3e2f1a0c9..."'
    """
    version_hash = get_model_version(*model_classes)
    return f'W/"{version_hash}"'


def generate_content_etag(content: str | bytes) -> str:
    """Generate a strong ETag from response content.

    Strong ETags indicate byte-for-byte equality. Use this when you need
    exact content matching.

    Parameters
    ----------
    content
        The response content to hash (string or bytes).

    Returns
    -------
    str
        Strong ETag in format: "<hash>"

    Examples
    --------
    >>> generate_content_etag("test content")
    '"9a0364b9e99bb480dd25e1f0284c8555"'
    >>> generate_content_etag(b"test content")
    '"9a0364b9e99bb480dd25e1f0284c8555"'
    """
    if isinstance(content, str):
        content = content.encode("utf-8")

    content_hash = hashlib.md5(content, usedforsecurity=False).hexdigest()
    return f'"{content_hash}"'


def parse_etag_header(header_value: str | None) -> list[str]:
    """Parse If-None-Match header into list of ETags.

    The If-None-Match header can contain:
    - A single ETag: "abc123"
    - Multiple ETags: "abc123", "def456"
    - Weak ETags: W/"abc123"
    - A wildcard: *

    Parameters
    ----------
    header_value
        The If-None-Match header value (or None).

    Returns
    -------
    list[str]
        List of parsed ETags. Empty list if header is None or empty.

    Examples
    --------
    >>> parse_etag_header('"abc123"')
    ['"abc123"']
    >>> parse_etag_header('"abc123", "def456"')
    ['"abc123"', '"def456"']
    >>> parse_etag_header('W/"abc123"')
    ['W/"abc123"']
    >>> parse_etag_header('*')
    ['*']
    >>> parse_etag_header(None)
    []
    """
    if not header_value:
        return []

    # Handle wildcard
    if header_value.strip() == "*":
        return ["*"]

    # Split by comma and strip whitespace
    etags = [etag.strip() for etag in header_value.split(",")]

    # Filter out empty strings
    return [etag for etag in etags if etag]
