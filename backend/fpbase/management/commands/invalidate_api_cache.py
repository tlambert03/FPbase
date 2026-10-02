from django.apps import apps
from django.core.management.base import BaseCommand

from fpbase import edge_cache
from fpbase.cache_utils import CACHED_MODELS, _invalidate
from fpbase.tasks import purge_edge_cache

# Run in the release phase, before the new web dyno starts: until it does, the old one
# can put pages and API responses back at the edge, so purge again once it is up.
DEPLOY_PURGE_DELAY = 120


class Command(BaseCommand):
    help = (
        "Drop every cached API response, on the server and at Cloudflare, and the cached "
        "pages at Cloudflare. Runs on every deploy; also for writes that send no signal: "
        "data migrations (run after `migrate`), raw SQL."
    )

    def handle(self, *args, **options):
        for label in sorted(CACHED_MODELS):
            _invalidate(apps.get_model(label))
        if edge_cache.is_enabled():
            purge_edge_cache.apply_async(kwargs={"pages": True}, countdown=DEPLOY_PURGE_DELAY)
            purged = "and Cloudflare purges are queued"
        else:
            purged = "(no CDN purge: not configured)"
        self.stdout.write(f"API caches invalidated {purged}")
