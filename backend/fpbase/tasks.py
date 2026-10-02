from __future__ import annotations

import logging

import requests
from celery import shared_task
from django.conf import settings

logger = logging.getLogger(__name__)

PURGE_PREFIXES = ("/api/", "/graphql/")
# the pages (besides the home page) in Cloudflare's "Cache public HTML" rule
PAGE_PURGE_PREFIXES = (
    "/protein/",
    "/reference/",
    "/spectra/",
    "/fret/",
    "/chart/",
    "/lineage/",
    "/microscopes/",
    "/organisms/",
    "/table/",
)


@shared_task(autoretry_for=(requests.RequestException,), retry_backoff=30, max_retries=2)
def purge_edge_cache(pages: bool = False) -> None:
    """Drop every cached API response at Cloudflare (see `fpbase.edge_cache`).

    With `pages`, also the HTML pages Cloudflare keeps for logged-out visitors.
    """
    prefixes = PURGE_PREFIXES + (PAGE_PURGE_PREFIXES if pages else ())
    host = settings.CANONICAL_URL.split("://", 1)[1]
    _purge({"prefixes": [host + prefix for prefix in prefixes]})
    if pages:
        # (a prefix of just the host would purge everything, static files included)
        _purge({"files": [settings.CANONICAL_URL + "/"]})


def _purge(payload: dict) -> None:
    response = requests.post(
        f"https://api.cloudflare.com/client/v4/zones/{settings.CLOUDFLARE_ZONE_ID}/purge_cache",
        json=payload,
        headers={"Authorization": f"Bearer {settings.CLOUDFLARE_PURGE_TOKEN}"},
        timeout=10,
    )
    if not response.ok:
        logger.error("Cloudflare purge failed: %s %s", response.status_code, response.text[:300])
