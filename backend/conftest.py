import json
import os
import random
import subprocess
from datetime import datetime
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from _pytest.tmpdir import TempPathFactory

REBUILD_ASSETS = "--rebuild-assets"


def pytest_addoption(parser):
    parser.addoption(
        REBUILD_ASSETS,
        action="store_true",
        default=False,
        help="skip building frontend",
    )


@pytest.fixture(scope="session")
def uses_frontend(request):
    from django.conf import settings

    stats_file = settings.WEBPACK_LOADER["DEFAULT"].get("STATS_FILE")
    if os.path.exists(stats_file):
        with open(stats_file, encoding="utf-8") as f:
            assets = json.load(f)
        assets_ready = assets.get("status") == "done" and assets.get("chunks")
    else:
        assets_ready = False

    if not assets_ready or request.config.getoption(REBUILD_ASSETS):
        subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)


@pytest.fixture(scope="session", autouse=True)
def _mock_blast_db(tmp_path_factory: "TempPathFactory"):
    from proteins.util import blast

    root = tmp_path_factory.mktemp("blastdb")

    blast.BLAST_DB, prev = str(root / "TEST_blastdb.fsa"), blast.BLAST_DB
    try:
        yield
    finally:
        blast.BLAST_DB = prev


@pytest.fixture()
def use_real_webpack_loader(monkeypatch):
    from webpack_loader import config, loaders, utils

    monkeypatch.setattr(
        utils,
        "get_loader",
        lambda config_name: loaders.WebpackLoader(config_name, config.load_config(config_name)),
    )


@pytest.fixture(autouse=True)
def mock_ncbi_api_calls(monkeypatch):
    """Mock NCBI API calls to avoid rate limiting during tests.

    This fixture automatically mocks external API calls to NCBI E-utilities
    to prevent HTTP 429 rate limit errors during test runs.
    """

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

    monkeypatch.setattr("proteins.extrest.entrez.doi_lookup", mock_doi_lookup)
    monkeypatch.setattr("proteins.extrest.entrez.get_organism_info", mock_get_organism_info)
