from __future__ import annotations

from django.apps import AppConfig


class FPbaseConfig(AppConfig):
    """Configuration for the fpbase app."""

    name = "fpbase"
    verbose_name = "FPbase"

    def ready(self):
        """Import signal handlers when Django starts."""
        from fpbase.cache_utils import _register_signal_handlers
        from proteins.models.microscope import get_cached_optical_configs
        from proteins.models.spectrum import get_cached_spectra_info

        _register_signal_handlers()

        # warm caches
        get_cached_spectra_info()
        get_cached_optical_configs()
