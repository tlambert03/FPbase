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
    from webpack_loader import config, loader, utils

    monkeypatch.setattr(
        utils,
        "get_loader",
        lambda config_name: loader.WebpackLoader(config_name, config.load_config(config_name)),
    )
