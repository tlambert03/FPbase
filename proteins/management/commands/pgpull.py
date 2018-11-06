from django.core.management.base import BaseCommand
from subprocess import run


class Command(BaseCommand):

    help = "wipe local database and pull from Heroku"

    def handle(self, *app_labels, **options):
        run(['dropdb', 'fpbase'])
        run(['heroku', 'pg:pull', 'DATABASE_URL', 'fpbase'])
        run(['python', 'manage.py', 'migrate'])
