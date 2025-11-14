"""Unified cache invalidation for model changes.

This module handles ALL cache invalidation across the application:

IMPORTANT: This is the ONLY place where post_save/post_delete signals
should be connected for cache invalidation purposes.
"""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING, Any

from django.apps import apps
from django.core.cache import cache
from django.db.models.signals import m2m_changed, post_delete, post_save
from django.dispatch import receiver
from django.utils import timezone

if TYPE_CHECKING:
    from django.db.models import Model


def _model_cache_key(model_class: type[Model]) -> str:
    return f"model_version:{model_class._meta.label}"


def get_model_version(*model_classes: type[Model]) -> str:
    """Get combined version hash for models."""
    versions = []
    for model_class in model_classes:
        cache_key = _model_cache_key(model_class)
        if (version := cache.get(cache_key)) is None:
            version = timezone.now().isoformat()
            if not cache.add(cache_key, version):
                # Another process set it first; get the value again
                version = cache.get(cache_key)
        versions.append(version)
    return hashlib.blake2b("".join(versions).encode(), digest_size=16).hexdigest()


def invalidate_model_version(model_class: type[Model]) -> None:
    """Bump the version for a model class."""
    cache.set(_model_cache_key(model_class), timezone.now().isoformat())


# Cache keys for JSON endpoints
SPECTRA_CACHE_KEY = "spectra_sluglist"
OPTICAL_CONFIG_CACHE_KEY = "optical_configs"


def _invalidate_spectra_cache() -> None:
    """Invalidate the spectra JSON cache."""
    cache.delete(SPECTRA_CACHE_KEY)


def _invalidate_optical_config_cache() -> None:
    """Invalidate the optical config JSON cache."""
    cache.delete(OPTICAL_CONFIG_CACHE_KEY)


SPECTRUM_OWNER_MODELS = {
    "proteins.Camera",
    "proteins.Dye",
    "proteins.Filter",
    "proteins.Light",
    "proteins.Protein",
    "proteins.Spectrum",
    "proteins.State",
}
OPTICAL_CONFIG_MODELS = {
    "proteins.Microscope",
    "proteins.OpticalConfig",
}


def _invalidate_on_change(sender: type[Model], **kwargs: Any) -> None:
    """Unified cache invalidation handler for model changes.

    This handler:
    1. Always invalidates model version (for ETags)
    2. Conditionally invalidates specific JSON caches based on model type
    """
    model_label = sender._meta.label

    # Always invalidate model version for ETag tracking
    invalidate_model_version(sender)

    # Invalidate spectra cache for models that affect spectra list
    if model_label in SPECTRUM_OWNER_MODELS:
        _invalidate_spectra_cache()

    # Invalidate optical config cache for models that affect optical configs
    if model_label in OPTICAL_CONFIG_MODELS:
        _invalidate_optical_config_cache()


def _register_signal_handlers():
    """Register signal handlers for model changes.

    IMPORTANT: This is the SINGLE place where these signals are connected.
    Must be called during app ready phase, not at module import time.
    """
    for model_label in SPECTRUM_OWNER_MODELS | OPTICAL_CONFIG_MODELS:
        model_class = apps.get_model(model_label)
        post_save.connect(
            _invalidate_on_change,
            sender=model_class,
            weak=False,
            dispatch_uid=f"{model_label}-cache-invalidate",
        )
        post_delete.connect(
            _invalidate_on_change,
            sender=model_class,
            weak=False,
            dispatch_uid=f"{model_label}-cache-invalidate-delete",
        )


@receiver(m2m_changed)
def invalidate_on_m2m_change(sender, instance, **kwargs):
    """Invalidate caches on many-to-many relationship changes."""
    invalidate_model_version(instance.__class__)
