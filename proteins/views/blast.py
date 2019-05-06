import json
import sys
import os
from django.core.files.storage import default_storage
from django.shortcuts import render
from django.http import JsonResponse
import tempfile
from Bio.Blast import NCBIXML
from subprocess import run


def serialize_alignment(alignment):
    out = alignment.__dict__.copy()
    out['hsps'] = [hsps.__dict__.copy() for hsps in alignment.hsps]
    return out


def serialize_record(record):
    out = record.__dict__.copy()
    out['database'] = 'FPbase'
    out['descriptions']
    out['alignments'] = [serialize_alignment(a) for a in record.alignments]
    # out['descriptions'] = [d.__dict__.copy() for d in record.descriptions]
    out.pop('descriptions')
    return out


def blastp(seq, binary='blastp', db='blastdb/FPbase_blastdb', max_hits=30, fmt=15, **kwargs):
    root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    binary = '{}/bin/{}_'.format(root, binary)
    binary += 'osx' if sys.platform == 'darwin' else 'nix'
    binary = os.path.abspath(binary)
    max_hits = kwargs.pop('max_target_seqs', max_hits)
    fmt = kwargs.pop('outfmt', fmt)
    with tempfile.NamedTemporaryFile(suffix='.fsa') as tmp:
        if not seq.startswith('>'):
            seq = '>query\n' + seq
        tmp.write(seq.encode())
        tmp.seek(0)
        cmd = [binary, '-query', tmp.name,
               '-outfmt', str(fmt),  # xml format
               '-db', '"' + default_storage.path(db) + '"',
               '-max_target_seqs', str(max_hits),
               '-max_hsps', '1'] # only show one alignment for each pair
        for key, value in kwargs.items():
            cmd.extend([f'-{key}', str(value)])
        with tempfile.NamedTemporaryFile(suffix='.txt') as outfile:
            cmd.extend(['-out', outfile.name])
            run(cmd)
            if fmt == 5:
                records = NCBIXML.parse(outfile.file)
                return [serialize_record(r) for r in records]
            elif fmt == 15:
                out = outfile.file.read().decode()
                return json.loads(out).get('BlastOutput2')
            else:
                return outfile.file.read().decode()


def blast_view(request):
    if request.is_ajax():
        seq = request.POST.get('query')
        binary = request.POST.get('binary', 'blastp')
        assert binary in ('blastx', 'blastp')
        if seq:
            return JsonResponse({
                'status': 200,
                'blastResult': blastp(seq, binary)
            })
        return JsonResponse({
            'status': 204
            })
    return render(request, 'proteins/blast.html')
