from __future__ import annotations

from django.apps import AppConfig


class FPbaseConfig(AppConfig):
    """Configuration for the fpbase app."""

    name = "fpbase"
    verbose_name = "FPbase"

    def ready(self):
        """Import signal handlers when Django starts."""
        import fpbase.cache_utils  # noqa: F401
