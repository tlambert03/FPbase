from __future__ import annotations

import logging

import requests
from celery import shared_task
from django.conf import settings

logger = logging.getLogger(__name__)

PURGE_PREFIXES = ("/api/", "/graphql/")


@shared_task(autoretry_for=(requests.RequestException,), retry_backoff=30, max_retries=2)
def purge_edge_cache() -> None:
    """Drop every cached API response at Cloudflare (see `fpbase.edge_cache`)."""
    host = settings.CANONICAL_URL.split("://", 1)[1]
    response = requests.post(
        f"https://api.cloudflare.com/client/v4/zones/{settings.CLOUDFLARE_ZONE_ID}/purge_cache",
        json={"prefixes": [host + prefix for prefix in PURGE_PREFIXES]},
        headers={"Authorization": f"Bearer {settings.CLOUDFLARE_PURGE_TOKEN}"},
        timeout=10,
    )
    if not response.ok:
        logger.error("Cloudflare purge failed: %s %s", response.status_code, response.text[:300])
