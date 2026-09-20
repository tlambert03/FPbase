import time

from django.core.cache import cache


def _protein_page_generation_key(slug: str) -> str:
    return f"protein-page-generation:{slug}"


def protein_page_key_prefix(slug: str) -> str:
    """`cache_page` key prefix for a protein's pages; changes when they are uncached."""
    return f"protein:{slug}:{cache.get(_protein_page_generation_key(slug), 0)}"


def uncache_protein_page(slug, request=None):
    # The page varies on cookies, so there is a cached copy per visitor and no way to
    # enumerate them: start a new generation rather than deleting (only) our own copy.
    cache.set(_protein_page_generation_key(slug), time.time_ns(), timeout=None)


def show_queries():
    import logging

    logger = logging.getLogger("django.db.backends")
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())


def is_ajax(request):
    # https://stackoverflow.com/a/70419609
    return request.headers.get("x-requested-with") == "XMLHttpRequest"
