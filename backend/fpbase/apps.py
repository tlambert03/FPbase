"""FPbase app configuration."""

from __future__ import annotations

from django.apps import AppConfig


class FPbaseConfig(AppConfig):
    """Configuration for the FPbase app."""

    name = "fpbase"
    verbose_name = "FPbase"
    default_auto_field = "django.db.models.BigAutoField"
