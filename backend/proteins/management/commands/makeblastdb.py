from django.core.management.base import BaseCommand

from proteins.util.blast import make_blastdb


class Command(BaseCommand):
    help = "Rebuild BLAST database"

    def handle(self, *app_labels, **options):
        make_blastdb()
