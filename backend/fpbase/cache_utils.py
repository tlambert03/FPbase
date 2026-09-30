"""Unified cache invalidation for model changes.

This module handles ALL cache invalidation across the application:

IMPORTANT: This is the ONLY place where post_save/post_delete signals
should be connected for cache invalidation purposes.
"""

from __future__ import annotations

import hashlib
import time
from functools import partial, wraps
from typing import TYPE_CHECKING, Any

from django.apps import apps
from django.core.cache import cache
from django.db import transaction
from django.db.models.signals import m2m_changed, post_delete, post_save
from django.dispatch import receiver
from django.utils import timezone
from django.utils.http import http_date
from django.views.decorators.cache import cache_page

if TYPE_CHECKING:
    from collections.abc import Callable

    from django.db.models import Model
    from django.http import HttpRequest, HttpResponse


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


DATA_VERSION_KEY = "data_version"
# Upper bound on how long a write that sends no signal (`queryset.update()`, raw SQL,
# a data migration) can go unnoticed by caches keyed on the data version.
DATA_CACHE_TTL = 60 * 60


def get_data_version() -> str:
    """Version of all the data served by the APIs; changes when any of it changes."""
    if (version := cache.get(DATA_VERSION_KEY)) is None:
        version = timezone.now().isoformat()
        if not cache.add(DATA_VERSION_KEY, version, None):
            version = cache.get(DATA_VERSION_KEY)
    return hashlib.blake2b(str(version).encode(), digest_size=8).hexdigest()


def cache_page_by_data_version(
    client_max_age: int = 600, public: bool = False
) -> Callable[[Callable[..., HttpResponse]], Callable[..., HttpResponse]]:
    """`cache_page`, with the data version in the key: a change to the data is a miss.

    Clients and the CDN know nothing of the data version, so they are given the
    shorter `client_max_age`, not the lifetime of the server-side cache entry.
    """
    cache_control = f"{'public, ' if public else ''}max-age={client_max_age}"

    def decorator(view: Callable[..., HttpResponse]) -> Callable[..., HttpResponse]:
        @wraps(view)
        def wrapper(request: HttpRequest, *args: Any, **kwargs: Any) -> HttpResponse:
            cached_view = cache_page(DATA_CACHE_TTL, key_prefix=get_data_version())(view)
            response = cached_view(request, *args, **kwargs)
            # cache_page stores a DRF response as it is rendered: do that first, so
            # that the headers of the stored copy are left as cache_page needs them
            if callable(render := getattr(response, "render", None)):
                render()
            if "max-age" in response.get("Cache-Control", ""):
                response["Cache-Control"] = cache_control
                response["Expires"] = http_date(time.time() + client_max_age)
                del response["Age"]
            return response

        return wrapper

    return decorator


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
    "proteins.DyeState",
    "proteins.Filter",
    "proteins.Light",
    "proteins.Protein",
    "proteins.Spectrum",
    "proteins.State",
}
OPTICAL_CONFIG_MODELS = {
    "proteins.Microscope",
    "proteins.OpticalConfig",
    "proteins.FilterPlacement",
}
# other models that appear in REST or GraphQL responses
API_MODELS = {
    "proteins.BleachMeasurement",
    "proteins.Excerpt",
    # a measurement rebuilds its state as a plain FluorState, not a State or DyeState
    "proteins.FluorescenceMeasurement",
    "proteins.FluorState",
    "proteins.Lineage",
    "proteins.OSERMeasurement",
    "proteins.StateTransition",
    "references.Author",
}
# models whose versions key the autocomplete search index (see proteins.search_index)
SEARCH_INDEX_MODELS = {
    "proteins.Dye",
    "proteins.DyeState",
    "proteins.Organism",
    "proteins.Protein",
    "proteins.Spectrum",
    "proteins.State",
    "references.Reference",
}


def _invalidate(sender: type[Model]) -> None:
    model_label = sender._meta.label

    invalidate_model_version(sender)
    if model_label in SPECTRUM_OWNER_MODELS:
        _invalidate_spectra_cache()
    if model_label in OPTICAL_CONFIG_MODELS:
        _invalidate_optical_config_cache()
    if model_label in CACHED_MODELS:
        cache.set(DATA_VERSION_KEY, timezone.now().isoformat(), None)


def _invalidate_on_change(sender: type[Model], **kwargs: Any) -> None:
    """Unified cache invalidation handler for model changes.

    This handler:
    1. Always invalidates model version (for ETags)
    2. Conditionally invalidates specific JSON caches based on model type
    3. Bumps the data version, if the model is served by the APIs
    """
    _invalidate(sender)
    if transaction.get_connection().in_atomic_block:
        # and again once the change is visible to other requests: one that read the
        # old rows in the meantime has cached them under the new version
        transaction.on_commit(partial(_invalidate, sender))


CACHED_MODELS = SPECTRUM_OWNER_MODELS | OPTICAL_CONFIG_MODELS | SEARCH_INDEX_MODELS | API_MODELS


def _register_signal_handlers():
    """Register signal handlers for model changes.

    IMPORTANT: This is the SINGLE place where these signals are connected.
    Must be called during app ready phase, not at module import time.
    """
    for model_label in CACHED_MODELS:
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
    _invalidate_on_change(sender=instance.__class__)
