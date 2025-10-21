"""NCBI API wrappers for Entrez, E-utilities, and related services.

All NCBI API calls should go through these functions to make mocking easier.
"""

from __future__ import annotations

import requests
from Bio import Entrez
from django.conf import settings

# Configure Entrez email/API key (used by all Bio.Entrez calls)
Entrez.email = getattr(settings, "NCBI_EMAIL", "talley.lambert+fpbase@gmail.org")
if api_key := getattr(settings, "NCBI_API_KEY", None):
    Entrez.api_key = api_key


# ============================================================================
# Bio.Entrez wrappers
# ============================================================================


def entrez_esearch(db: str, term: str, **kwargs):
    """Search NCBI database using Entrez.esearch.

    Parameters
    ----------
    db : str
        Database to search (e.g., 'pubmed', 'taxonomy', 'ipg', 'protein')
    term : str
        Search term
    **kwargs
        Additional arguments to pass to Entrez.esearch

    Returns
    -------
    Bio.Entrez result handle

    """
    return Entrez.esearch(db=db, term=term, **kwargs)


def entrez_esummary(db: str, id: str | int, retmode: str = "xml", **kwargs):
    """Get document summaries using Entrez.esummary.

    Parameters
    ----------
    db : str
        Database to query (e.g., 'pubmed', 'taxonomy', 'ipg')
    id : str | int
        Record ID(s) to retrieve
    retmode : str
        Return mode ('xml' or 'json')
    **kwargs
        Additional arguments to pass to Entrez.esummary

    Returns
    -------
    Bio.Entrez result handle

    """
    return Entrez.esummary(db=db, id=id, retmode=retmode, **kwargs)


def entrez_efetch(db: str, id: str | int, rettype: str | None = None, retmode: str = "text", **kwargs):
    """Fetch full records using Entrez.efetch.

    Parameters
    ----------
    db : str
        Database to fetch from (e.g., 'protein', 'nuccore')
    id : str | int
        Record ID(s) to fetch
    rettype : str | None
        Return type (e.g., 'fasta', 'gb')
    retmode : str
        Return mode (e.g., 'text', 'xml')
    **kwargs
        Additional arguments to pass to Entrez.efetch

    Returns
    -------
    Bio.Entrez result handle

    """
    return Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode, **kwargs)


def entrez_espell(db: str, term: str, **kwargs):
    """Get spelling suggestions using Entrez.espell.

    Parameters
    ----------
    db : str
        Database to check spelling against
    term : str
        Term to check
    **kwargs
        Additional arguments to pass to Entrez.espell

    Returns
    -------
    Bio.Entrez result handle

    """
    return Entrez.espell(db=db, term=term, **kwargs)


def entrez_read(handle):
    """Read and parse an Entrez result handle.

    Parameters
    ----------
    handle
        Entrez result handle to read

    Returns
    -------
    Parsed data structure (usually dict or list)

    """
    return Entrez.read(handle)


# ============================================================================
# Direct E-utilities REST API wrappers (via requests)
# ============================================================================


def eutils_efetch(
    db: str,
    id: str | int,
    rettype: str = "fasta",
    retmode: str = "text",
    **params,
) -> requests.Response:
    """Fetch records using E-utilities REST API directly.

    Parameters
    ----------
    db : str
        Database to fetch from (e.g., 'protein', 'nuccore', 'pmc')
    id : str | int
        Record ID to fetch
    rettype : str
        Return type (e.g., 'fasta', 'gb')
    retmode : str
        Return mode (e.g., 'text', 'xml')
    **params
        Additional query parameters

    Returns
    -------
    requests.Response
        HTTP response from E-utilities

    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params_dict = {
        "db": db,
        "id": id,
        "rettype": rettype,
        "retmode": retmode,
        **params,
    }
    return requests.get(base_url, params=params_dict)


def pmc_id_converter(id: str, tool: str = "FPbase", email: str | None = None) -> dict:
    """Convert between PubMed, PMC, and DOI identifiers.

    Uses NCBI's PMC ID Converter API.

    Parameters
    ----------
    id : str
        Identifier to convert (PMID, PMCID, or DOI)
    tool : str
        Tool name for API tracking
    email : str | None
        Email for API tracking

    Returns
    -------
    dict
        Conversion result with 'pmid', 'pmcid', and 'doi' keys

    """
    email = email or Entrez.email
    url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool={tool}&email={email}&ids={id}&format=json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()
