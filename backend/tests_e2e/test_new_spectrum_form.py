from __future__ import annotations

import re
from typing import TYPE_CHECKING

from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import ProteinFactory
from proteins.models import Spectrum

if TYPE_CHECKING:
    from collections.abc import Callable

    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


def _tomselect_enter(selector: str, text: str, page: Page) -> None:
    """Helper to select an option in a Tom-Select widget by typing and selecting."""
    # Tom-Select creates an input field with class "ts-control"
    # Find the Tom-Select input associated with the original select element
    ts_input = page.locator(f"{selector} + .ts-wrapper .ts-control input")
    ts_input.click()
    ts_input.fill(text)
    # Wait for dropdown options to appear
    option = page.locator(".ts-dropdown .option").first
    expect(option).to_be_visible()
    option.click()


def test_spectrum_submission_preview_manual_data(
    auth_page: Page, live_server: LiveServer, assert_snapshot: Callable
) -> None:
    """Test spectrum submission form with manual data preview."""
    protein = ProteinFactory.create()
    protein.default_state.ex_spectrum.delete()

    url = f"{live_server.url}/spectra-modern/submit/"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)
    # Wait for form to be fully initialized
    expect(auth_page.locator("#spectrum-form")).to_be_attached()

    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)
    auth_page.locator("#id_confirmation").check()

    # Tom-Select autocomplete for protein
    _tomselect_enter("#id_owner_state", protein.name, auth_page)

    # Switch to manual data tab and enter data
    auth_page.locator("#manual-tab").click()
    data_field = auth_page.locator("#id_spectral_data")
    expect(data_field).to_be_visible()
    data_field.fill("[[500,0.1],[505,0.5],[510,0.8],[515,0.6],[520,0.3]]")

    # Visual snapshot: form filled but before preview
    assert_snapshot(auth_page)

    # Submit for preview
    auth_page.locator('button[type="submit"]').click()

    # Wait for preview section to appear (AJAX request)
    preview_section = auth_page.locator("#spectrum-preview-section")
    expect(preview_section).to_be_visible(timeout=10000)

    svg = auth_page.locator("#spectrum-preview-chart svg")
    expect(svg).to_be_visible()
    expect(svg.locator("[id^='FillBetweenPolyCollection']")).to_have_count(1)

    # Visual snapshot: preview chart displayed
    if not hasattr(assert_snapshot, "NOOP"):
        auth_page.wait_for_load_state("networkidle")
        assert_snapshot(auth_page)

    # submit it!
    auth_page.get_by_text("Submit Spectrum").click()
    expect(auth_page).to_have_url(f"{live_server.url}{reverse('spectra_modern:submitted')}")


def test_spectrum_submission_tab_switching(
    auth_page: Page, live_server: LiveServer, assert_snapshot: Callable
) -> None:
    """Test tab switching behavior in spectrum submission form."""
    # Create a protein so owner_state field has options
    protein = ProteinFactory.create(name="TabSwitchTestProtein")

    url = f"{live_server.url}/spectra-modern/submit/"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Fill out basic fields
    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)
    auth_page.locator("#id_confirmation").check()

    # Wait for Tom-Select to initialize on the protein field (it initializes when category changes to 'p')
    auth_page.wait_for_selector("#id_owner_state + .ts-wrapper", state="visible")
    _tomselect_enter("#id_owner_state", protein.name, auth_page)

    # Test tab switching
    file_tab = auth_page.locator("#file-tab")
    manual_tab = auth_page.locator("#manual-tab")

    # Verify file tab is active by default
    expect(file_tab).to_have_class(re.compile("active"))

    # Switch to manual tab
    manual_tab.click()
    expect(manual_tab).to_have_class(re.compile("active"))

    # Enter manual data
    data_field = auth_page.locator("#id_spectral_data")
    expect(data_field).to_be_visible()
    expect(data_field).to_be_enabled()
    data_field.fill("[[400,0.1],[401,0.2],[402,0.3],[403,0.5],[404,0.8],[405,1.0]]")

    # Visual snapshot: manual tab active with data
    expect(auth_page.get_by_text("File Upload")).not_to_be_visible()
    assert_snapshot(auth_page)

    # Switch back to file tab
    file_tab.click()
    expect(file_tab).to_have_class(re.compile("active"))

    # Verify submit button is present
    submit_btn = auth_page.locator('button[type="submit"]')
    expect(submit_btn).to_be_visible()

    # Visual snapshot: file tab active
    expect(auth_page.get_by_text("File Upload")).to_be_visible()
    assert_snapshot(auth_page)
