"""Model version tracking for cache invalidation.

This module provides utilities for tracking when Django models change,
enabling efficient cache invalidation and ETag generation.
"""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

from django.core.cache import cache
from django.db.models.signals import m2m_changed, post_delete, post_save
from django.dispatch import receiver
from django.utils import timezone

if TYPE_CHECKING:
    from django.db.models import Model


def get_model_version(*model_classes: type[Model]) -> str:
    """Get combined version hash for one or more models.

    Each model's version is tracked in cache. When a model instance is saved
    or deleted, the version is automatically bumped via signals.

    Parameters
    ----------
    *model_classes
        One or more Django model classes to get versions for.

    Returns
    -------
    str
        MD5 hash of combined model versions. If only one model, returns
        hash of that model's version. If multiple, returns hash of all
        versions concatenated.

    Examples
    --------
    >>> get_model_version(Spectrum)
    'a3c2f1e8b9d4...'
    >>> get_model_version(Spectrum, Protein)
    'b4d3e2f1a0c9...'
    """
    versions = []
    for model_class in model_classes:
        # Use app_label.ModelName for cache key
        cache_key = f"model_version:{model_class._meta.label}"
        version = cache.get(cache_key)

        # Initialize version if not in cache
        if version is None:
            version = timezone.now().isoformat()
            cache.set(cache_key, version)

        versions.append(version)

    # Combine all versions and hash
    combined = "".join(versions)
    return hashlib.md5(combined.encode()).hexdigest()


def invalidate_model_version(model_class: type[Model]) -> None:
    """Bump the version for a model class.

    This is automatically called by signal handlers when instances are
    saved or deleted. You generally don't need to call this manually.

    Parameters
    ----------
    model_class
        Django model class to invalidate.

    Examples
    --------
    >>> invalidate_model_version(Spectrum)
    """
    cache_key = f"model_version:{model_class._meta.label}"
    new_version = timezone.now().isoformat()
    cache.set(cache_key, new_version)


# Signal handlers for automatic version invalidation
# Register these for models that need version tracking


@receiver([post_save, post_delete], sender="proteins.Spectrum")
def invalidate_spectrum_version(sender, **kwargs):
    """Invalidate Spectrum version on save/delete."""
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Protein")
def invalidate_protein_version(sender, **kwargs):
    """Invalidate Protein version on save/delete."""
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.State")
def invalidate_state_version(sender, **kwargs):
    """Invalidate State version on save/delete."""
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.OpticalConfig")
def invalidate_optical_config_version(sender, **kwargs):
    """Invalidate OpticalConfig version on save/delete."""
    invalidate_model_version(sender)


# Handle many-to-many changes (e.g., protein.states.add())
@receiver(m2m_changed)
def invalidate_on_m2m_change(sender, instance, **kwargs):
    """Invalidate version when m2m relationships change."""
    invalidate_model_version(instance.__class__)
