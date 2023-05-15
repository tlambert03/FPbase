import sys
from subprocess import run

from django.core.management.base import BaseCommand

from proteins.models import Protein


class Command(BaseCommand):
    help = "Cache Forster distance calculations across database"

    def handle(self, *app_labels, **options):
        seqs = Protein.objects.filter(parent_organism=6100)
        seqs = seqs.exclude(name__in=["GCaMP2", "cEGFP", "mKalama1", "mAmetrine"])
        with open("fasta.fasta", "w") as handle:
            handle.write(seqs.fasta().read())
        if sys.platform == "darwin":
            binary = "bin/muscle_osx"
        else:
            binary = "bin/muscle_nix"
        cmd = [binary]
        # faster
        cmd += ["-maxiters", "2", "-diags", "-quiet", "-clw", "-out", "test.clw"]
        cmd += ["-cluster", "neighborjoining", "-tree2", "tree.phy"]
        p = run(cmd, input=seqs.fasta().read(), encoding="ascii")
        print(p.stdout)
