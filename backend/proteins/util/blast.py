import contextlib
import json
import os
import shutil
import sys
import tempfile
from functools import cache
from pathlib import Path
from shutil import copyfileobj
from subprocess import run

from ..models import Protein

ROOT = Path(__file__).parent.parent.parent
BIN_DIR = ROOT / "bin"
BLAST_DB = "blastdb/FPbase_blastdb.fsa"
BIN_SUFFIX = "osx" if sys.platform == "darwin" else "nix"


@cache
def _get_binary(binary: str) -> str:
    """Get the path to a binary, or raise an error if it's not found"""
    local_version = BIN_DIR / f"{binary}_{BIN_SUFFIX}"
    if local_version.exists():
        return str(local_version)

    fallback = shutil.which(binary)
    if fallback:
        return fallback
    raise FileNotFoundError(f"{binary} program not found")


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


def write_fasta(fpath):
    """Writes all FPsequences to fasta file in default storage location"""

    os.makedirs(os.path.dirname(fpath), exist_ok=True)

    with open(fpath, "w") as fd:
        with contextlib.suppress(Exception):
            # for some reason, first write usually throws an exception
            fd.write("")
        fasta = Protein.objects.all().fasta()
        fasta.seek(0)
        copyfileobj(fasta, fd)
        return fd.name


def make_blastdb(fpath: str | None = None):
    fasta_name = write_fasta(fpath or BLAST_DB)
    cmd = [
        _get_binary("makeblastdb"),
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


def blast(seq, binary="blastp", db: str | None = None, max_hits=30, fmt=15, **kwargs):
    db = db or BLAST_DB

    assert binary in ("blastp", "blastx"), "Unrecognized blast binary"
    binary = _get_binary(binary)
    if not os.path.isfile(db) or len(os.listdir(os.path.dirname(db))) <= 5:
        make_blastdb(db)

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
            return _extracted_from_blast_(cmd, outfile, fmt)


# TODO Rename this here and in `blast`
def _extracted_from_blast_(cmd, outfile, fmt):
    cmd.extend(["-out", outfile.name])
    run(cmd)
    if fmt == 5:
        from Bio.Blast import NCBIXML

        records = NCBIXML.parse(outfile.file)
        return [serialize_record(r) for r in records]
    if fmt != 15:
        return outfile.file.read().decode()
    out = outfile.file.read().decode()
    return json.loads(out).get("BlastOutput2")
