from django.core.management.base import BaseCommand
from subprocess import call


class Command(BaseCommand):

    help = "Start completely fresh."

    def handle(self, *app_labels, **options):
        call(['dropdb', 'fpbase'])
        call(['createdb', 'fpbase'])
        call(['python', 'manage.py', 'migrate'])
        call(['python', 'manage.py', 'loaddata', 'users'])
        #call(['python', 'manage.py', 'loaddata', 'proteins', '-i'])
        #call(['python', 'manage.py', 'createinitialrevisions'])
