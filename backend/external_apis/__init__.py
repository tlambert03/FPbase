"""Centralized external API wrappers.

This package provides thin wrappers around all external API calls used in FPbase.
By centralizing these calls, we make it easier to:
1. Mock external services during testing
2. Track and manage rate limits
3. Add retry logic and error handling consistently
4. Switch implementations if needed

All external API calls should go through these wrappers rather than calling
external libraries directly.
"""

from __future__ import annotations

__all__ = [
    "ncbi",
    "references",
    "sequences",
]
