from django.core.management.base import BaseCommand
from subprocess import run, PIPE
from proteins.models import Protein


class Command(BaseCommand):

    help = "Cache Forster distance calculations across database"

    def handle(self, *app_labels, **options):
        seqs = Protein.objects.filter(parent_organism=6100).fasta()
        p = run('bin/muscle_nix', stdout=PIPE, input=seqs.read(), encoding='ascii')
        print(p.stdout)
