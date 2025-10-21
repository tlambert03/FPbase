"""External APIs for bibliographic references (Crossref, etc.).

All reference/publication API calls should go through these functions.
"""

from __future__ import annotations

from habanero import Crossref


def crossref_works(doi: str, mailto: str = "talley.lambert+fpbase@gmail.org") -> dict:
    """Fetch publication metadata from Crossref API.

    Parameters
    ----------
    doi : str
        Digital Object Identifier to look up
    mailto : str
        Email address for Crossref API (enables polite pool)

    Returns
    -------
    dict
        Crossref work metadata

    """
    cr = Crossref(mailto=mailto)
    return cr.works(ids=doi)
