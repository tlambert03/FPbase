"""App configuration for modern spectrum submission."""

from django.apps import AppConfig


class SpectraModernConfig(AppConfig):
    """Configuration for the spectra_modern app."""

    name = "spectra_modern"
    verbose_name = "Modern Spectrum Submission"
    default_auto_field = "django.db.models.BigAutoField"
