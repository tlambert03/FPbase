"""Model version tracking for cache invalidation."""

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
    """Get combined version hash for models."""
    versions = []
    for model_class in model_classes:
        cache_key = f"model_version:{model_class._meta.label}"
        version = cache.get(cache_key)
        if version is None:
            version = timezone.now().isoformat()
            cache.set(cache_key, version)
        versions.append(version)
    return hashlib.md5("".join(versions).encode()).hexdigest()


def invalidate_model_version(model_class: type[Model]) -> None:
    """Bump the version for a model class."""
    cache_key = f"model_version:{model_class._meta.label}"
    cache.set(cache_key, timezone.now().isoformat())


@receiver([post_save, post_delete], sender="proteins.Spectrum")
def invalidate_spectrum_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Protein")
def invalidate_protein_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.State")
def invalidate_state_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.OpticalConfig")
def invalidate_optical_config_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Microscope")
def invalidate_microscope_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Dye")
def invalidate_dye_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Camera")
def invalidate_camera_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Light")
def invalidate_light_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver([post_save, post_delete], sender="proteins.Filter")
def invalidate_filter_version(sender, **kwargs):
    invalidate_model_version(sender)


@receiver(m2m_changed)
def invalidate_on_m2m_change(sender, instance, **kwargs):
    invalidate_model_version(instance.__class__)
