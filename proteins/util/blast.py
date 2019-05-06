import sys
import os
from ..models import Protein
from .helpers import remember_cwd
from django.core.files.storage import default_storage
from subprocess import run


def write_fasta(fpath='blastdb/FPbase_blastdb.fsa'):
    """Writes all FPsequences to fasta file in default storage location"""
    from shutil import copyfileobj

    with default_storage.open(fpath, 'w') as fd:
        try:
            # for some reason, first write usually throws an exception
            fd.write('')
        except Exception:
            pass
        fasta = Protein.objects.all().fasta()
        fasta.seek(0)
        copyfileobj(fasta, fd)
        return fd.name


def make_blastdb(dirname='blastdb', name='FPbase_blastdb'):
    root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    binary = root + '/bin/makeblastdb_' + 'osx' if sys.platform == 'darwin' else 'nix'
    binary = os.path.abspath(binary)
    with remember_cwd():
        os.chdir(default_storage.path(dirname))
        write_fasta(default_storage.path(os.path.join(dirname, name + '.fsa')))
        cmd = [binary, '-in', name + '.fsa', '-parse_seqids', '-blastdb_version',
               '5', '-title', 'FPbase Sequence Database', '-dbtype', 'prot',
               '-out', name]
        run(cmd)
