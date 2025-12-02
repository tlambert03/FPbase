"""Run fpbasepy's test suite against our live server with seeded test data."""

import os
import subprocess

import pytest

from proteins.factories import (
    CameraFactory,
    DyeStateFactory,
    FilterFactory,
    LightFactory,
    MicroscopeFactory,
    OpticalConfigFactory,
    ProteinFactory,
    create_egfp,
)
from proteins.models import Filter

FPBASEPY_REPO = "https://github.com/tlambert03/fpbasepy"


@pytest.fixture
def seed_fpbasepy_data(db):
    """Seed database with data expected by fpbasepy tests."""
    # Proteins expected by tests
    create_egfp()
    ProteinFactory(
        name="mScarlet-I",
        default_state__name="default",
        default_state__ex_max=569,
        default_state__em_max=593,
        default_state__ext_coeff=104000,
        default_state__qy=0.54,
    )
    ProteinFactory(
        name="mEos3.2",
        default_state__name="default",
        default_state__ex_max=507,
        default_state__em_max=516,
        default_state__ext_coeff=63400,
        default_state__qy=0.84,
    )
    # Additional proteins for pdb test
    ProteinFactory(name="Clover1.5")
    ProteinFactory(name="6C")
    ProteinFactory(name="dClover2 A206K")

    # Dye expected by tests
    DyeStateFactory(dye__name="Alexa Fluor 488", name="default")

    # Filters expected by tests
    FilterFactory(name="Chroma ET525/50m", subtype=Filter.BPM)
    FilterFactory(name="Semrock FF01-520/35", subtype=Filter.BP)

    # Camera and light source
    CameraFactory(name="Andor Zyla 5.5")
    LightFactory(name="Lumencor Celesta UV")

    # Microscope with optical configs
    scope = MicroscopeFactory(name="Example Simple Widefield", id="wKqWbgApvguSNDSRZNSfpN")
    ex_filter = FilterFactory(name="Chroma ET470/40x", subtype=Filter.BPX)
    bs_filter = FilterFactory(name="Chroma T495lpxr", subtype=Filter.LP)
    em_filter = FilterFactory(name="Chroma ET525/50m-2", subtype=Filter.BPM)
    oc = OpticalConfigFactory(name="Widefield Green", microscope=scope)
    oc.filters.add(ex_filter, through_defaults={"path": "EX"})
    oc.filters.add(bs_filter, through_defaults={"path": "BS"})
    oc.filters.add(em_filter, through_defaults={"path": "EM"})


@pytest.mark.django_db(transaction=True)
def test_fpbasepy_upstream_suite(live_server, tmp_path, seed_fpbasepy_data):
    """Clone and run fpbasepy's test suite against our GraphQL API."""
    repo = tmp_path / "fpbasepy"

    # Clone and install
    subprocess.run(["git", "clone", "--depth=1", FPBASEPY_REPO, str(repo)], check=True)
    env = {k: v for k, v in os.environ.items() if k != "UV_FROZEN"}
    subprocess.run(["uv", "sync", "--extra", "test", "-q"], cwd=repo, env=env, check=True)

    # Patch URL to point at live server
    (repo / "tests" / "conftest.py").write_text(f"""\
import fpbase._fetch
fpbase._fetch.FPbaseClient._FPbaseClient__instance = fpbase._fetch.FPbaseClient(
    base_url="{live_server.url}/graphql/"
)
""")

    # Run tests (exclude test_get_missing_protein - it tests fpbasepy's fuzzy matching)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("VIRTUAL_ENV", "DJANGO"))}
    result = subprocess.run(
        ["uv", "run", "pytest", "tests", "-v", "-s", "-k", "not test_get_missing"],
        cwd=repo,
        env=env,
    )
    assert result.returncode == 0, "fpbasepy test suite failed"
