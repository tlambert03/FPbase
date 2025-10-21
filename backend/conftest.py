import json
import os
import subprocess
from io import BytesIO
from typing import TYPE_CHECKING

import pytest
from django.conf import settings

from proteins.util import blast

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
def mock_external_apis(monkeypatch):
    """Mock all external API calls to avoid rate limiting and network dependencies during tests.

    This fixture mocks the centralized external_apis module, making it much easier
    to maintain than mocking individual Bio.Entrez, habanero, and requests calls.
    """
    # Mock NCBI API wrapper functions
    def mock_entrez_esearch(db, term, **kwargs):
        """Mock NCBI esearch to return a test record handle."""
        xml = b"""<?xml version="1.0" encoding="UTF-8"?>
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
        return BytesIO(xml)

    def mock_entrez_esummary(db, id, retmode="xml", **kwargs):
        """Mock NCBI esummary to return test data."""
        if db == "pubmed":
            xml = b"""<?xml version="1.0" encoding="UTF-8"?>
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
        elif db == "taxonomy":
            xml = b"""<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN"
            "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
            <eSummaryResult>
                <DocSum>
                    <Id>9606</Id>
                    <Item Name="Rank" Type="String">species</Item>
                    <Item Name="Division" Type="String">hydrozoans</Item>
                    <Item Name="ScientificName" Type="String">Aequorea victoria</Item>
                    <Item Name="CommonName" Type="String">jellyfish</Item>
                    <Item Name="Species" Type="String">victoria</Item>
                    <Item Name="Genus" Type="String">Aequorea</Item>
                </DocSum>
            </eSummaryResult>"""
        else:
            xml = b"""<?xml version="1.0" encoding="UTF-8"?>
            <eSummaryResult><DocSum></DocSum></eSummaryResult>"""
        return BytesIO(xml)

    def mock_entrez_efetch(db, id, rettype=None, retmode="text", **kwargs):
        """Mock NCBI efetch to return test sequence data."""
        return BytesIO(b">test\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS")

    def mock_entrez_espell(db, term, **kwargs):
        """Mock NCBI espell to return the same term."""
        xml = b"""<?xml version="1.0" encoding="UTF-8"?>
        <eSpellResult>
            <Query>""" + term.encode() + b"""</Query>
            <CorrectedQuery></CorrectedQuery>
        </eSpellResult>"""
        return BytesIO(xml)

    def mock_entrez_read(handle):
        """Mock Entrez.read - just pass through to Bio.Entrez.read."""
        from Bio import Entrez

        return Entrez.read(handle)

    def mock_pmc_id_converter(id, tool="FPbase", email=None):
        """Mock PMC ID converter to return test data."""
        return {
            "records": [
                {
                    "pmid": "12345678",
                    "pmcid": "PMC1234567",
                    "doi": "10.1234/test.doi",
                }
            ]
        }

    # Mock Crossref API wrapper
    def mock_crossref_works(doi, mailto="talley.lambert+fpbase@gmail.org"):
        """Mock Crossref works API to return test publication data."""
        return {
            "message": {
                "DOI": doi,
                "title": [f"Test Article for {doi}"],
                "container-title": ["Test Journal of Science"],
                "page": "123-456",
                "volume": "42",
                "issue": "1",
                "published-print": {"date-parts": [[2024, 1, 1]]},
                "author": [
                    {"family": "Smith", "given": "John A."},
                    {"family": "Doe", "given": "Jane B."},
                ],
            }
        }

    # Mock sequence fetcher APIs
    def mock_fetch_genbank_fasta(accession):
        """Mock GenBank FASTA fetch."""
        return ">test\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    def mock_fetch_uniprot_fasta(accession):
        """Mock UniProt FASTA fetch."""
        return ">test\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    def mock_fetch_pdb_sequence(accession):
        """Mock PDB sequence fetch."""
        return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"

    def mock_fetch_uniprot_xml(uniprot_id):
        """Mock UniProt XML fetch."""
        return """<?xml version="1.0" encoding="UTF-8"?>
        <uniprot>
            <entry>
                <accession>TEST123</accession>
                <sequence>MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS</sequence>
            </entry>
        </uniprot>"""

    def mock_map_uniprot_ids(ids, from_db, to_db):
        """Mock UniProt ID mapping."""
        return "\n".join([f"{id}\tTEST{i}" for i, id in enumerate(ids)])

    # Apply all mocks to the centralized external_apis module
    monkeypatch.setattr("external_apis.ncbi.entrez_esearch", mock_entrez_esearch)
    monkeypatch.setattr("external_apis.ncbi.entrez_esummary", mock_entrez_esummary)
    monkeypatch.setattr("external_apis.ncbi.entrez_efetch", mock_entrez_efetch)
    monkeypatch.setattr("external_apis.ncbi.entrez_espell", mock_entrez_espell)
    monkeypatch.setattr("external_apis.ncbi.entrez_read", mock_entrez_read)
    monkeypatch.setattr("external_apis.ncbi.pmc_id_converter", mock_pmc_id_converter)
    monkeypatch.setattr("external_apis.references.crossref_works", mock_crossref_works)
    monkeypatch.setattr("external_apis.sequences.fetch_genbank_fasta", mock_fetch_genbank_fasta)
    monkeypatch.setattr("external_apis.sequences.fetch_uniprot_fasta", mock_fetch_uniprot_fasta)
    monkeypatch.setattr("external_apis.sequences.fetch_pdb_sequence", mock_fetch_pdb_sequence)
    monkeypatch.setattr("external_apis.sequences.fetch_uniprot_xml", mock_fetch_uniprot_xml)
    monkeypatch.setattr("external_apis.sequences.map_uniprot_ids", mock_map_uniprot_ids)

    # Also mock references.models imports (used directly in some places)
    monkeypatch.setattr("references.models.name_to_initials", lambda x: x[:2].upper())

    # Mock Organism.save() to skip NCBI API calls
    from proteins.models.organism import Organism

    def mock_organism_save(self, *args, **kwargs):
        # Set default values if not already set, skipping API call
        if not self.scientific_name:
            self.scientific_name = "Aequorea victoria"
        if not self.division:
            self.division = "hydrozoans"
        if not self.common_name:
            self.common_name = "jellyfish"
        if not self.species:
            self.species = "victoria"
        if not self.genus:
            self.genus = "Aequorea"
        if not self.rank:
            self.rank = "species"
        # Call Django's Model.save() directly, skipping Organism's custom save
        super(Organism, self).save(*args, **kwargs)

    monkeypatch.setattr(Organism, "save", mock_organism_save)
