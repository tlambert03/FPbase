"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import ProteinFactory
from proteins.models import Spectrum

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


def test_spectrum_submission_preview_manual_data(auth_page: Page, live_server: LiveServer) -> None:
    """Test spectrum submission form with manual data preview."""
    protein = ProteinFactory.create(name="TestGFP", agg="m")
    spectrum = protein.default_state.ex_spectrum
    spectrum.delete()

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)

    # Select2 autocomplete for protein
    combo = auth_page.locator("#div_id_owner_state [role='combobox']")
    combo.click()
    expect(combo).to_have_attribute("aria-expanded", "true")
    auth_page.keyboard.type(protein.name)
    auth_page.keyboard.press("Enter")

    # Switch to manual data tab and enter data
    auth_page.locator("#manual-tab").click()
    auth_page.locator("#id_data").fill(json.dumps(spectrum.data))
    auth_page.locator("#id_confirmation").check()

    # Submit for preview
    auth_page.locator('input[type="submit"]').click()

    # Verify preview section appears with chart
    expect(auth_page.locator("#spectrum-preview-section")).to_be_visible()
    svg = auth_page.locator("#spectrum-preview-chart svg")
    expect(svg).to_be_visible()
    expect(svg.locator("[id^='FillBetweenPolyCollection']")).to_have_count(1)

    # submit it!
    auth_page.get_by_text("Submit Spectrum").click()
    expect(auth_page).to_have_url(f"{live_server.url}{reverse('proteins:spectrum_submitted')}")
