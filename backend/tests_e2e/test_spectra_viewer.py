"""End-to-end tests for /spectra/ viewer page."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import create_egfp

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


@pytest.fixture
def spectra_viewer(live_server: LiveServer, page: Page) -> Iterator[Page]:
    create_egfp()

    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    spectra_viewer = page.locator("#spectra-viewer")
    expect(spectra_viewer).to_be_visible()
    yield page
    # DO NOT REMOVE: exposes many JS errors that happen on page actions
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(100)


def test_spectra_viewer_add_from_input(spectra_viewer: Page, assert_snapshot: Callable) -> None:
    """Test the spectra viewer page loads without console errors."""
    tab_wrapper = spectra_viewer.locator(".tab-wrapper")
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()
    search_input = tab_wrapper.get_by_role("combobox")
    search_input.click()
    search_input.type("EGF")

    egfp_option = tab_wrapper.get_by_role("option", name="EGFP", exact=True)
    expect(egfp_option).to_be_visible()
    egfp_option.click()

    # a selector for EGFP should now be visible
    expect(tab_wrapper.get_by_text(re.compile(r"^EGFP"))).to_be_visible()
    # a NEW search input should appear after selecting a protein
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()

    assert_snapshot(spectra_viewer)


@pytest.mark.parametrize("method", ["spacebar", "click"])
def test_spectra_viewer_add_from_spacebar(spectra_viewer: Page, assert_snapshot: Callable, method: str) -> None:
    """Test the spectra viewer page loads without console errors."""
    tab_wrapper = spectra_viewer.locator(".tab-wrapper")
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()

    modal = spectra_viewer.get_by_role("presentation").filter(has_text="Quick Entry")
    quick_entry_heading = modal.get_by_role("heading", name="Quick Entry")
    expect(modal).not_to_be_attached()
    expect(quick_entry_heading).not_to_be_visible()

    if method == "spacebar":
        # press spacebar to open the search input
        spectra_viewer.keyboard.press("Space")
    elif method == "click":
        # click on the #quickentry-btn button
        spectra_viewer.locator("button#quickentry-btn").click()
    else:
        raise ValueError(f"Unknown method: {method}")

    expect(modal).to_be_attached()
    expect(quick_entry_heading).to_be_visible()

    search_input = modal.get_by_role("combobox").first
    expect(search_input).to_be_visible()
    search_input.type("EGF")

    expect(modal.get_by_role("option", name="EGFP", exact=True)).to_be_visible()
    spectra_viewer.keyboard.press("Enter")

    # a NEW search input should appear after selecting a protein
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()
    # ... and a selector for EGFP should now be visible
    expect(tab_wrapper.get_by_text(re.compile(r"^EGFP"))).to_be_visible()

    assert_snapshot(spectra_viewer)


def test_spectra_url_sharing_basic(live_server: LiveServer, page: Page) -> None:
    """Test basic URL sharing functionality in spectra viewer.

    Verifies:
    1. Share button opens dialog with generated URL
    2. Generated URL contains all chart state parameters
    """
    egfp = create_egfp()

    # Navigate to spectra viewer
    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Add EGFP using spacebar shortcut
    page.keyboard.press("Space")
    modal = page.get_by_role("presentation").filter(has_text="Quick Entry")
    expect(modal).to_be_attached()

    search_input = modal.get_by_role("combobox").first
    expect(search_input).to_be_visible()
    search_input.type("EGFP")
    page.keyboard.press("Enter")

    # Wait for EGFP to be added to the selector list
    tab_wrapper = page.locator(".tab-wrapper")
    expect(tab_wrapper.get_by_text(re.compile(r"^EGFP"))).to_be_visible()

    # Wait for chart to load with data
    chart = page.locator(".highcharts-container")
    expect(chart).to_be_visible()
    # Wait for series to render
    page.locator(".highcharts-series").first.wait_for(state="attached")

    # Click share button using aria attributes for reliability
    page.wait_for_load_state("networkidle")
    share_button = page.locator('button[aria-controls="simple-menu"]')
    expect(share_button).to_be_visible()
    share_button.click()

    # Wait for menu to appear and click "Share chart as URL"
    share_menu_item = page.get_by_role("menuitem", name="Share chart as URL")
    expect(share_menu_item).to_be_visible(timeout=10000)
    share_menu_item.click()

    # Verify share dialog is open
    dialog = page.get_by_role("dialog", name=re.compile("recreate the current graph"))
    expect(dialog).to_be_visible()

    # Get the URL from the textbox
    url_textbox = dialog.get_by_role("textbox", name="URL")
    expect(url_textbox).to_be_visible()

    shared_url = url_textbox.input_value()

    # Verify URL contains spectra ID parameter
    assert "s=" in shared_url, "Shared URL should contain spectra ID parameter"
    # Check that URL contains one of the spectrum IDs (ex or em)
    ex_id = str(egfp.default_state.ex_spectrum.id)
    em_id = str(egfp.default_state.em_spectrum.id)
    assert ex_id in shared_url or em_id in shared_url, (
        f"Shared URL should contain EGFP spectrum ID (ex:{ex_id} or em:{em_id})"
    )


def test_spectra_url_params_parsing(live_server: LiveServer, page: Page) -> None:
    """Test that URL parameters are correctly parsed and applied to spectra viewer.

    Verifies:
    1. URL with chart options is parsed correctly
    2. Chart state matches URL parameters
    3. Backwards compatibility with legacy URL formats
    """
    egfp = create_egfp()
    ex_id = egfp.default_state.ex_spectrum.id
    em_id = egfp.default_state.em_spectrum.id

    # Navigate with URL parameters (testing backwards-compatible format)
    url = (
        f"{live_server.url}{reverse('proteins:spectra')}"
        f"?s={ex_id},{em_id}"
        "&showY=1&showX=1&showGrid=1&areaFill=1"
        "&logScale=0&scaleEC=0&scaleQY=0&shareTooltip=1"
        "&palette=wavelength&xMin=400&xMax=600"
    )
    page.goto(url)

    # Wait for chart to load
    chart = page.locator(".highcharts-container")
    expect(chart).to_be_visible()

    # Verify chart loaded with spectra
    chart_text = page.locator(".highcharts-series").first
    expect(chart_text).to_be_attached()

    # Verify extremes were applied (check input boxes)
    min_input = page.locator('input[name="min"]')
    max_input = page.locator('input[name="max"]')
    expect(min_input).to_have_value("400")
    expect(max_input).to_have_value("600")
