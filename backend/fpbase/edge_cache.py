"""API responses cached at the CDN (Cloudflare), purged when the data changes.

Off until `CLOUDFLARE_ZONE_ID` and `CLOUDFLARE_PURGE_TOKEN` are set: the CDN is only
told to keep responses for `EDGE_MAX_AGE` when they will also be purged.
"""

from __future__ import annotations

from django.conf import settings
from django.core.cache import cache

from fpbase.tasks import purge_edge_cache

# (a change purges sooner; see also `manage.py invalidate_api_cache`, run on every deploy)
EDGE_MAX_AGE = 24 * 60 * 60
BROWSER_MAX_AGE = 60
# changes committed within this many seconds share one purge (Cloudflare allows 5/min)
PURGE_DELAY = 10
PURGE_PENDING_KEY = "edge_purge_pending"


def is_enabled() -> bool:
    return bool(settings.CLOUDFLARE_ZONE_ID and settings.CLOUDFLARE_PURGE_TOKEN)


def cache_control(public: bool = False, client_max_age: int = 600) -> str:
    """The Cache-Control header for a response that is valid until the data changes."""
    if is_enabled():
        return f"public, max-age={BROWSER_MAX_AGE}, s-maxage={EDGE_MAX_AGE}"
    return f"{'public, ' if public else ''}max-age={client_max_age}"


def schedule_purge() -> None:
    """Purge the API paths at the CDN in `PURGE_DELAY` seconds, unless one is pending."""
    if not is_enabled():
        return
    # (`add` is None, not False, when the cache is unreachable: purge anyway)
    if cache.add(PURGE_PENDING_KEY, 1, PURGE_DELAY) is not False:
        purge_edge_cache.apply_async(countdown=PURGE_DELAY)
