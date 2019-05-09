import sys
from ..models import Protein
from subprocess import run
import os
import tempfile
import json


def serialize_alignment(alignment):
    out = alignment.__dict__.copy()
    out["hsps"] = [hsps.__dict__.copy() for hsps in alignment.hsps]
    return out


def serialize_record(record):
    out = record.__dict__.copy()
    out["database"] = "FPbase"
    out["descriptions"]
    out["alignments"] = [serialize_alignment(a) for a in record.alignments]
    # out['descriptions'] = [d.__dict__.copy() for d in record.descriptions]
    out.pop("descriptions")
    return out


def write_fasta(fpath="blastdb/FPbase_blastdb.fsa"):
    """Writes all FPsequences to fasta file in default storage location"""
    from shutil import copyfileobj

    os.makedirs(os.path.dirname(fpath), exist_ok=True)

    with open(fpath, "w") as fd:
        try:
            # for some reason, first write usually throws an exception
            fd.write("")
        except Exception:
            pass
        fasta = Protein.objects.all().fasta()
        fasta.seek(0)
        copyfileobj(fasta, fd)
        return fd.name


def make_blastdb(fpath="blastdb/FPbase_blastdb.fsa"):
    binary = "bin/makeblastdb_" + ("osx" if sys.platform == "darwin" else "nix")
    fasta_name = write_fasta(fpath)
    cmd = [
        binary,
        "-in",
        fasta_name,
        "-parse_seqids",
        "-blastdb_version",
        "5",
        "-title",
        "FPbase Sequence Database",
        "-dbtype",
        "prot",
    ]
    run(cmd)


def blast(
    seq, binary="blastp", db="blastdb/FPbase_blastdb.fsa", max_hits=30, fmt=15, **kwargs
):
    assert binary in ("blastp", "blastx"), "Unrecognized blast binary"
    if not (os.path.isfile(db) and (len(os.listdir(os.path.dirname(db))) > 5)):
        make_blastdb(db)
    binary = f"bin/{binary}_" + ("osx" if sys.platform == "darwin" else "nix")
    max_hits = kwargs.pop("max_target_seqs", max_hits)
    fmt = kwargs.pop("outfmt", fmt)
    with tempfile.NamedTemporaryFile(suffix=".fsa") as tmp:
        if not seq.startswith(">"):
            seq = ">query\n" + seq
        tmp.write(seq.encode())
        tmp.seek(0)
        cmd = [
            binary,
            "-query",
            tmp.name,
            "-outfmt",
            str(fmt),  # xml format
            "-db",
            db,
            "-max_target_seqs",
            str(max_hits),
            "-max_hsps",
            "1",
        ]  # only show one alignment for each pair
        for key, value in kwargs.items():
            cmd.extend([f"-{key}", str(value)])
        with tempfile.NamedTemporaryFile(suffix=".txt") as outfile:
            cmd.extend(["-out", outfile.name])
            run(cmd)
            if fmt == 5:
                from Bio.Blast import NCBIXML

                records = NCBIXML.parse(outfile.file)
                return [serialize_record(r) for r in records]
            if fmt == 15:
                out = outfile.file.read().decode()
                return json.loads(out).get("BlastOutput2")
            else:
                return outfile.file.read().decode()

