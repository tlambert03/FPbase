from django.core.cache import cache
from proteins.util.helpers import forster_list
from django.core.management.base import BaseCommand


class Command(BaseCommand):

    help = "Cache Forster distance calculations across database"

    def handle(self, *app_labels, **options):
        L = forster_list()
        cache.set("forster_list", L, 60 * 60 * 24)
