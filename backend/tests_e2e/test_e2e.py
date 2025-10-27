"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import ProteinFactory

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


@pytest.mark.usefixtures("assert_no_console_errors")
def test_main_page_loads_with_assets(live_server: LiveServer, page: Page) -> None:
    """Test that the main page loads without errors and CSS/JS assets are applied.

    Verifies:
    1. Page loads without console errors
    2. Navigation elements are visible and accessible
    3. CSS styling is applied (via semantic checks, not computed styles)
    4. Interactive elements are present

    This is a foundational test demonstrating Playwright best practices:
    - Uses semantic locators (get_by_role, get_by_placeholder)
    - Auto-waiting via expect() assertions
    - Console error detection via fixture
    - No explicit waits or time.sleep()
    """
    # Navigate to main page
    page.goto(live_server.url)

    # Verify navigation is visible and accessible
    # Using get_by_role for semantic locator (WCAG compliant)
    navbar = page.get_by_role("navigation")
    expect(navbar).to_be_visible()

    # Verify logo link is present (semantic check - looks for link to home)
    # This checks both presence and that it's a functioning link
    home_link = page.get_by_role("link", name="FPbase")
    expect(home_link).to_be_visible()

    # Verify search functionality is present
    # Uses placeholder text for semantic identification
    search_input = page.get_by_placeholder("Search")
    expect(search_input).to_be_visible()
    expect(search_input).to_be_enabled()

    # Verify main content area exists
    # Rather than checking computed styles, verify semantic structure
    main_content = page.locator("main, #content, .main-content")
    expect(main_content.first).to_be_visible()


@pytest.mark.usefixtures("assert_no_console_errors")
def test_spectra_viewer_loads(live_server: LiveServer, page: Page) -> None:
    """Test the spectra viewer page loads without console errors."""
    ProteinFactory(name="TestGFP", agg="m")

    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    spectra_viewer = page.locator("#spectra-viewer")
    expect(spectra_viewer).to_be_attached()
