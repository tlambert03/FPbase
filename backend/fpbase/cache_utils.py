"""Model version tracking for cache invalidation."""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

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
            cache.set(cache_key, version)
        versions.append(version)
    return hashlib.blake2b("".join(versions).encode(), digest_size=16).hexdigest()


def invalidate_model_version(model_class: type[Model]) -> None:
    """Bump the version for a model class."""
    cache.set(_model_cache_key(model_class), timezone.now().isoformat())


def _invalidate_on_change(sender, **kwargs):
    invalidate_model_version(sender)


# Register signal handlers for model changes
for model_label in [
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
    model_class = apps.get_model(model_label)
    post_save.connect(
        _invalidate_on_change,
        sender=model_class,
        weak=False,
        dispatch_uid=f"{model_label}-post_save-invalidate",
    )
    post_delete.connect(
        _invalidate_on_change,
        sender=model_class,
        weak=False,
        dispatch_uid=f"{model_label}-post_delete-invalidate",
    )


@receiver(m2m_changed)
def invalidate_on_m2m_change(sender, instance, **kwargs):
    invalidate_model_version(instance.__class__)
