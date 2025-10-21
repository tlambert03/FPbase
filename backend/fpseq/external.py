from __future__ import annotations

from external_apis import sequences


def genbank_seq(accession: str) -> str | None:
    """Retrieve protein sequence from GenBank using NCBI E-utilities.

    Parameters
    ----------
    accession : str
        GenBank accession ID

    Returns
    -------
    str | None
        Protein sequence, or None if not found
    """
    try:
        fasta_text = sequences.fetch_genbank_fasta(accession)
        # Parse FASTA: skip header line, join sequence lines
        return "".join(fasta_text.splitlines()[1:])
    except Exception:
        return None


def uniprot_seq(accession: str) -> str | None:
    """Retrieve protein sequence from UniProt.

    Parameters
    ----------
    accession : str
        UniProt accession ID

    Returns
    -------
    str | None
        Protein sequence, or None if not found
    """
    try:
        fasta_text = sequences.fetch_uniprot_fasta(accession)
        # Parse FASTA: skip header line, join sequence lines
        return "".join(fasta_text.splitlines()[1:])
    except Exception:
        return None


def pdb_seq(accession: str) -> str | None:
    """Retrieve protein sequence from PDB.

    Parameters
    ----------
    accession : str
        PDB accession ID

    Returns
    -------
    str | None
        Protein sequence, or None if not found
    """
    try:
        return sequences.fetch_pdb_sequence(accession)
    except Exception:
        return None
