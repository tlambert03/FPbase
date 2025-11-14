from __future__ import annotations

from django.apps import AppConfig


class FPbaseConfig(AppConfig):
    """Configuration for the fpbase app."""

    name = "fpbase"
    verbose_name = "FPbase"

    def ready(self):
        """Import signal handlers when Django starts."""
        from fpbase.cache_utils import _register_signal_handlers

        _register_signal_handlers()
