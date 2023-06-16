import logging
import subprocess
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from _pytest.tmpdir import TempPathFactory


def pytest_addoption(parser):
    parser.addoption(
        "--skip-build",
        action="store_true",
        default=False,
        help="skip building frontend",
    )


@pytest.fixture(scope="session")
def uses_frontend(request):
    if not request.config.getoption("--skip-build"):
        logging.info("Building frontend...")
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
