import json
import os
import subprocess
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
    from datetime import datetime

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
            "pmid": "12345678",
        }

    def mock_get_pmid_info(pmid):
        """Mock get_pmid_info to return realistic test data."""
        return {
            "doi": "10.1234/test.doi",
            "title": f"Test Article for PMID {pmid}",
            "journal": "Test Journal",
            "pages": "123-456",
            "volume": "42",
            "issue": "1",
            "year": "2024",
            "authors": [{"family": "Smith", "given": "John A."}],
            "date": datetime(2024, 1, 1).date(),
        }

    def mock_doi2pmid(doi):
        """Mock doi2pmid to return a test PMID."""
        return "12345678"

    def mock_pmid2doi(pmid):
        """Mock pmid2doi to return a test DOI."""
        return "10.1234/test.doi"

    def mock_entrez_esearch(db, term, **kwargs):
        """Mock Entrez.esearch to return a test record."""
        from io import StringIO

        # Return XML that Bio.Entrez.read can parse
        xml = """<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN"
        "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
        <eSearchResult>
            <Count>1</Count>
            <RetMax>1</RetMax>
            <RetStart>0</RetStart>
            <IdList>
                <Id>12345678</Id>
            </IdList>
        </eSearchResult>"""
        return StringIO(xml)

    def mock_entrez_esummary(db, id, retmode="xml", **kwargs):
        """Mock Entrez.esummary to return test data."""
        from io import StringIO

        if db == "pubmed":
            xml = """<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN"
            "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
            <eSummaryResult>
                <DocSum>
                    <Id>12345678</Id>
                    <Item Name="PubDate" Type="Date">2024 Jan 01</Item>
                    <Item Name="EPubDate" Type="Date">2024 Jan 01</Item>
                    <Item Name="Source" Type="String">Test Journal</Item>
                    <Item Name="Title" Type="String">Test Article Title</Item>
                    <Item Name="Volume" Type="String">42</Item>
                    <Item Name="Issue" Type="String">1</Item>
                    <Item Name="Pages" Type="String">123-456</Item>
                    <Item Name="DOI" Type="String">10.1234/test.doi</Item>
                    <Item Name="AuthorList" Type="List">
                        <Item Name="Author" Type="String">Smith JA</Item>
                        <Item Name="Author" Type="String">Doe JB</Item>
                    </Item>
                </DocSum>
            </eSummaryResult>"""
        else:
            xml = """<?xml version="1.0" encoding="UTF-8"?>
            <eSummaryResult><DocSum></DocSum></eSummaryResult>"""
        return StringIO(xml)

    # Apply the mocks
    monkeypatch.setattr("references.helpers.doi_lookup", mock_doi_lookup)
    monkeypatch.setattr("references.helpers.get_pmid_info", mock_get_pmid_info)
    monkeypatch.setattr("references.helpers.doi2pmid", mock_doi2pmid)
    monkeypatch.setattr("references.helpers.pmid2doi", mock_pmid2doi)

    def mock_genbank_seq(accession):
        """Mock genbank_seq to return a test sequence."""
        return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    def mock_uniprot_seq(accession):
        """Mock uniprot_seq to return a test sequence."""
        return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    def mock_pdb_seq(accession):
        """Mock pdb_seq to return a test sequence."""
        return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    # Apply the mocks
    monkeypatch.setattr("references.helpers.doi_lookup", mock_doi_lookup)
    monkeypatch.setattr("references.helpers.get_pmid_info", mock_get_pmid_info)
    monkeypatch.setattr("references.helpers.doi2pmid", mock_doi2pmid)
    monkeypatch.setattr("references.helpers.pmid2doi", mock_pmid2doi)

    # Mock external sequence fetchers
    monkeypatch.setattr("fpseq.external.genbank_seq", mock_genbank_seq)
    monkeypatch.setattr("fpseq.external.uniprot_seq", mock_uniprot_seq)
    monkeypatch.setattr("fpseq.external.pdb_seq", mock_pdb_seq)

    # Mock Bio.Entrez functions
    try:
        from Bio import Entrez

        monkeypatch.setattr(Entrez, "esearch", mock_entrez_esearch)
        monkeypatch.setattr(Entrez, "esummary", mock_entrez_esummary)
    except ImportError:
        pass  # Bio not imported yet
