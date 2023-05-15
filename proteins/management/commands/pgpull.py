from subprocess import run

from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = "wipe local database and pull from Heroku"

    def handle(self, *app_labels, **options):
        run(["dropdb", "fpbase"])
        run(["heroku", "pg:pull", "DATABASE_URL", "fpbase", "-a", "fpbase"])
        run(["python", "manage.py", "migrate"])
