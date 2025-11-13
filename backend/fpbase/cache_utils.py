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


def _cache_key(model_class: type[Model]) -> str:
    return f"model_version:{model_class._meta.label}"


def get_model_version(*model_classes: type[Model]) -> str:
    """Get combined version hash for models."""
    versions = []
    for model_class in model_classes:
        cache_key = _cache_key(model_class)
        if (version := cache.get(cache_key)) is None:
            cache.set(cache_key, timezone.now().isoformat())
        versions.append(version)
    return hashlib.blake2b("".join(versions).encode(), digest_size=16).hexdigest()


def invalidate_model_version(model_class: type[Model]) -> None:
    """Bump the version for a model class."""
    cache.set(_cache_key(model_class), timezone.now().isoformat())


# Register signal handlers for model changes
for model_name in [
    "proteins.Spectrum",
    "proteins.Protein",
    "proteins.State",
    "proteins.OpticalConfig",
    "proteins.Microscope",
    "proteins.Dye",
    "proteins.Camera",
    "proteins.Light",
    "proteins.Filter",
]:
    post_save.connect(lambda sender, **_: invalidate_model_version(sender), sender=model_name)
    post_delete.connect(lambda sender, **_: invalidate_model_version(sender), sender=model_name)


@receiver(m2m_changed)
def invalidate_on_m2m_change(sender, instance, **kwargs):
    invalidate_model_version(instance.__class__)
