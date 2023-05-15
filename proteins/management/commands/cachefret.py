from django.core.cache import cache
from django.core.management.base import BaseCommand

from proteins.util.helpers import forster_list


class Command(BaseCommand):
    help = "Cache Forster distance calculations across database"

    def handle(self, *app_labels, **options):
        L = forster_list()
        cache.set("forster_list", L, 60 * 60 * 24)
