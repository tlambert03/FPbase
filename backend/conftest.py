import random
from datetime import datetime
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from _pytest.tmpdir import TempPathFactory


@pytest.fixture(scope="session", autouse=True)
def _mock_blast_db(tmp_path_factory: "TempPathFactory"):
    from proteins.util import blast

    root = tmp_path_factory.mktemp("blastdb")

    blast.BLAST_DB, prev = str(root / "TEST_blastdb.fsa"), blast.BLAST_DB
    try:
        yield
    finally:
        blast.BLAST_DB = prev


@pytest.fixture(scope="session", autouse=True)
def mock_ncbi_api_calls():
    """Mock NCBI API calls to avoid rate limiting during tests.

    This fixture automatically mocks external API calls to NCBI E-utilities
    to prevent HTTP 429 rate limit errors during test runs.

    Session-scoped to ensure mocks are active during Django TestCase.setUpTestData().
    """
    import unittest.mock

    def mock_doi_lookup(doi):
        """Mock doi_lookup to return realistic test data without API calls."""
        return {
            "title": f"Test Article for {doi}",
            "journal": "Test Journal of Science",
            "pages": "123-456",
            "volume": "42",
            "issue": "1",
            "year": 2024,
            "date": datetime(2024, 1, 1).date(),
            "authors": [
                {"family": "Smith", "given": "John A."},
                {"family": "Doe", "given": "Jane B."},
            ],
            # random 8 numbers
            "pmid": random.randint(10_000_000, 99_999_999),
        }

    def mock_get_organism_info(organism_id: str | int):
        return {
            "scientific_name": "ScientificName",
            "division": "Division",
            "common_name": "CommonName",
            "species": "Species",
            "genus": "Genus",
            "rank": "Rank",
        }

    # Use unittest.mock.patch for session-scoped mocking
    patcher1 = unittest.mock.patch(
        "proteins.extrest.entrez.doi_lookup", side_effect=mock_doi_lookup
    )
    patcher2 = unittest.mock.patch(
        "proteins.extrest.entrez.get_organism_info", side_effect=mock_get_organism_info
    )

    patcher1.start()
    patcher2.start()

    yield

    patcher1.stop()
    patcher2.stop()
