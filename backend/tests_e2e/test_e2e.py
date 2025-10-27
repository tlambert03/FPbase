"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


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
