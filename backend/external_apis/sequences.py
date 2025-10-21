"""External APIs for protein/nucleotide sequence retrieval.

All sequence fetching from external databases should go through these functions.
"""

from __future__ import annotations

import requests


def fetch_genbank_fasta(accession: str) -> str:
    """Fetch protein sequence from GenBank in FASTA format.

    Uses NCBI E-utilities REST API directly.

    Parameters
    ----------
    accession : str
        GenBank accession number

    Returns
    -------
    str
        FASTA-formatted sequence

    Raises
    ------
    requests.HTTPError
        If the request fails

    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={accession}&rettype=fasta&retmode=text"
    response = requests.get(url)
    response.raise_for_status()
    return response.text


def fetch_uniprot_fasta(accession: str) -> str:
    """Fetch protein sequence from UniProt in FASTA format.

    Parameters
    ----------
    accession : str
        UniProt accession number

    Returns
    -------
    str
        FASTA-formatted sequence

    Raises
    ------
    requests.HTTPError
        If the request fails

    """
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    response = requests.get(url)
    response.raise_for_status()
    return response.text


def fetch_pdb_sequence(accession: str) -> str:
    """Fetch protein sequence from Protein Data Bank via EBI API.

    Parameters
    ----------
    accession : str
        PDB accession code

    Returns
    -------
    str
        Protein sequence

    Raises
    ------
    requests.HTTPError
        If the request fails
    ValueError
        If the PDB entry has no sequence data

    """
    url = f"http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{accession}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    # Extract the first chain's sequence
    if data.get(accession):
        for molecule in data[accession]:
            if "sequence" in molecule:
                return molecule["sequence"]

    raise ValueError(f"No sequence found for PDB entry {accession}")


def fetch_uniprot_xml(uniprot_id: str) -> str:
    """Fetch full UniProt record in XML format.

    Parameters
    ----------
    uniprot_id : str
        UniProt ID or accession

    Returns
    -------
    str
        XML-formatted UniProt record

    Raises
    ------
    requests.HTTPError
        If the request fails

    """
    url = f"http://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    response.raise_for_status()
    return response.text


def map_uniprot_ids(ids: list[str], from_db: str, to_db: str) -> str:
    """Map between UniProt and other database identifiers.

    Parameters
    ----------
    ids : list[str]
        List of identifiers to map
    from_db : str
        Source database code (e.g., 'ACC' for UniProt accession)
    to_db : str
        Target database code (e.g., 'EMBL' for GenBank)

    Returns
    -------
    str
        Mapping results in tab-separated format

    Raises
    ------
    requests.HTTPError
        If the request fails

    """
    url = "http://www.uniprot.org/uploadlists/"
    params = {
        "from": from_db,
        "to": to_db,
        "format": "tab",
        "query": " ".join(ids),
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.text
