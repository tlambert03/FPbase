from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import django.conf
import pytest

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer

# Required for pytest-playwright fixtures to work with Django ORM in test setup.
# Playwright runs fixtures in an async event loop, triggering Django's async-unsafe
# protections. This is safe for tests because:
# 1) test fixtures run sequentially
# 2) live_server runs Django normally in a separate thread
# 3) no concurrent ORM access
# See: https://docs.djangoproject.com/en/5.2/topics/async/#envvar-DJANGO_ALLOW_ASYNC_UNSAFE
#      https://github.com/microsoft/playwright-pytest/issues/29
os.environ.setdefault("DJANGO_ALLOW_ASYNC_UNSAFE", "true")


@pytest.fixture(scope="module", autouse=True)
def _build_frontend_assets() -> None:
    """Build webpack assets once per test session if needed.

    Checks if existing webpack stats represent a static build.
    Rebuilds if missing, stale, or --rebuild-assets flag is set.
    """
    stats_file = Path(django.conf.settings.WEBPACK_LOADER["DEFAULT"]["STATS_FILE"])
    if stats_file.is_file():
        assets = json.loads(stats_file.read_bytes())
        if (
            assets.get("status") == "done"
            and assets.get("chunks")
            and ("localhost" not in assets.get("publicPath", ""))
        ):
            return

    print("Building frontend assets for e2e tests...")
    subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)


def test_main_page(live_server: LiveServer, page: Page) -> None:
    """Test the main page loads and CSS assets are applied."""
    page.goto(live_server.url)

    # Check that CSS was loaded by verifying computed styles
    # FPbase uses Bootstrap, so check for common Bootstrap styling
    body = page.locator("body")
    font_family = body.evaluate("el => window.getComputedStyle(el).fontFamily")

    # Bootstrap sets a specific font stack - if it's just default, CSS didn't load
    assert font_family and font_family != "Times New Roman", f"CSS not loaded properly, got font: {font_family}"

    # Verify the custom FPbase background gradient is applied (shows CSS loaded)
    body_bg = body.evaluate("el => window.getComputedStyle(el).backgroundImage")
    assert "gradient" in body_bg.lower(), f"Background gradient not applied: {body_bg}"

    # Verify styled elements are present
    navbar = page.locator("nav.navbar")
    assert navbar.count() > 0, "Navbar not found"

    # Check for the FPbase logo
    logo = page.locator('a[href="/"] img, a[href="/"] svg')
    assert logo.count() > 0, "FPbase logo not found"

    # Verify search box is styled (has proper Bootstrap classes or custom styling)
    search_input = page.locator('input[type="search"], input[placeholder*="Search"]')
    if search_input.count() > 0:
        padding = search_input.first.evaluate("el => window.getComputedStyle(el).padding")
        # If CSS loaded, padding should be set (not default "0px")
        assert padding != "0px", f"Search input not styled, padding: {padding}"
