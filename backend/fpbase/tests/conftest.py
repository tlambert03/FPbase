import logging
import subprocess

import pytest


def pytest_addoption(parser):
    parser.addoption("--skip-build", action="store_true", default=False, help="skip building frontend")


@pytest.fixture(scope="session")
def uses_frontend(request):
    if not request.config.getoption("--skip-build"):
        logging.info("Building frontend...")
        subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)


@pytest.fixture()
def use_real_webpack_loader(monkeypatch):
    from django.conf import settings

    monkeypatch.delitem(settings.WEBPACK_LOADER["DEFAULT"], "LOADER_CLASS")
