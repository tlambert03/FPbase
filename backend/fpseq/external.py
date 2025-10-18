from __future__ import annotations

import json
import os

import requests


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
    api_key = os.getenv("NCBI_API_KEY", "")
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=protein&id={accession}&rettype=fasta&retmode=text"
    )
    if api_key:
        url += f"&api_key={api_key}"

    response = requests.get(url)
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
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
    response = requests.get(f"https://www.uniprot.org/uniprot/{accession}.fasta")
    if response:
        return "".join([x.decode() for x in response.content.splitlines()[1:]])
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
    response = requests.get(f"http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{accession}")
    if response:
        j = json.loads(response.text)
        try:
            j = j[accession.lower()]
        except ValueError:
            j = j[accession]
        return j[0]["sequence"]
    return None
