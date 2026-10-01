from django.apps import apps
from django.core.management.base import BaseCommand

from fpbase import edge_cache
from fpbase.cache_utils import CACHED_MODELS, _invalidate


class Command(BaseCommand):
    help = (
        "Drop every cached API response, on the server and at Cloudflare. For writes that "
        "send no signal: data migrations (run after `migrate`), raw SQL."
    )

    def handle(self, *args, **options):
        for label in sorted(CACHED_MODELS):
            _invalidate(apps.get_model(label))
        purged = (
            "and a Cloudflare purge is queued"
            if edge_cache.is_enabled()
            else "(no CDN purge: not configured)"
        )
        self.stdout.write(f"API caches invalidated {purged}")
