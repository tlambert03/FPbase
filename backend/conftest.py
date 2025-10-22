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
    needs_rebuild = True

    if os.path.exists(stats_file):
        with open(stats_file, encoding="utf-8") as f:
            assets = json.load(f)

        # Check if assets are from a production build (not dev server)
        # Dev builds have publicPath pointing to localhost:8080
        is_production_build = assets.get("status") == "done" and assets.get("chunks")
        if is_production_build and assets.get("publicPath"):
            # If publicPath contains localhost, it's from dev server - need rebuild
            is_production_build = "localhost" not in assets.get("publicPath")

        needs_rebuild = not is_production_build

    if needs_rebuild or request.config.getoption(REBUILD_ASSETS):
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
    patcher1 = unittest.mock.patch("proteins.extrest.entrez.doi_lookup", side_effect=mock_doi_lookup)
    patcher2 = unittest.mock.patch("proteins.extrest.entrez.get_organism_info", side_effect=mock_get_organism_info)

    patcher1.start()
    patcher2.start()

    yield

    patcher1.stop()
    patcher2.stop()
