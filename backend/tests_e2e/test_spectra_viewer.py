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


def _add_egfp_to_viewer(page: Page) -> None:
    tab_wrapper = page.locator(".tab-wrapper")
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()
    search_input = tab_wrapper.get_by_role("combobox")
    search_input.click()
    search_input.type("EGFP")
    egfp_option = tab_wrapper.get_by_role("option", name="EGFP", exact=True)
    expect(egfp_option).to_be_visible()
    egfp_option.click()
    expect(tab_wrapper.get_by_text(re.compile(r"^EGFP"))).to_be_visible()


def test_spectra_viewer_add_from_input(spectra_viewer: Page, assert_snapshot: Callable) -> None:
    """Test the spectra viewer page loads without console errors."""
    _add_egfp_to_viewer(spectra_viewer)
    # a NEW search input should appear after selecting a protein
    tab_wrapper = spectra_viewer.locator(".tab-wrapper")
    expect(tab_wrapper.get_by_text("Type to search...")).to_be_visible()
    if not hasattr(assert_snapshot, "NOOP"):
        # Wait for chart to fully render
        spectra_viewer.locator(".highcharts-series").first.wait_for(state="attached")
        spectra_viewer.wait_for_load_state("networkidle")
        # Mask legend as it has minor rendering variations
        assert_snapshot(spectra_viewer, mask_elements=[".highcharts-legend"])


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

    if not hasattr(assert_snapshot, "NOOP"):
        # Wait for chart to fully render
        spectra_viewer.locator(".highcharts-series").first.wait_for(state="attached")
        spectra_viewer.wait_for_load_state("networkidle")
        # Mask legend as it has minor rendering variations
        assert_snapshot(spectra_viewer, mask_elements=[".highcharts-legend"])


def test_spectra_url_sharing_basic(live_server: LiveServer, page: Page) -> None:
    """Test basic URL sharing functionality in spectra viewer.

    Verifies:
    1. Share button opens dialog with generated URL
    2. Generated URL contains all chart state parameters
    """
    create_egfp()

    # Navigate to spectra viewer
    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Add EGFP to the viewer
    _add_egfp_to_viewer(page)

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

    # Verify URL contains spectra ID parameter with at least one ID
    assert "s=" in shared_url, "Shared URL should contain spectra ID parameter"
    # Extract spectrum IDs from URL (format: ?s=123,456 or ?s=123)
    match = re.search(r"[?&]s=([0-9,]+)", shared_url)
    assert match, f"Could not find spectrum IDs in URL: {shared_url}"
    spectrum_ids = match.group(1).split(",")
    assert len(spectrum_ids) > 0, "Shared URL should contain at least one spectrum ID"
    # Verify IDs are numeric
    assert all(sid.isdigit() for sid in spectrum_ids), f"Invalid spectrum IDs: {spectrum_ids}"


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


def test_qy_ec_scaling_invertibility(spectra_viewer: Page) -> None:
    """Test that QY and EC scaling transformations are invertible.

    Regression test for bug where toggling QY/EC scaling multiple times
    would compound transformations instead of applying them idempotently.
    """

    spectra_viewer.keyboard.press("Space")
    modal = spectra_viewer.get_by_role("presentation").filter(has_text="Quick Entry")
    search_input = modal.get_by_role("combobox").first
    search_input.type("EGFP")
    spectra_viewer.keyboard.press("Enter")

    # Wait for EGFP to be added
    expect(spectra_viewer.locator(".tab-wrapper").get_by_text(re.compile(r"^EGFP"))).to_be_visible()

    def _test_scaling_toggle_indempotent(page: Page, series: int, key: str) -> None:
        _ser = page.locator(f"g.highcharts-series.highcharts-series-{series}")
        line = _ser.locator("path.highcharts-tracker-line")
        # Wait for line to be visible AND have data (not empty/loading)
        expect(line).to_be_visible()
        expect(line).to_have_attribute("d", re.compile(r"^M .+"))

        # Get initial path data
        d_0 = line.get_attribute("d")
        assert d_0 is not None, "Initial path data should not be None"

        # Press key once and wait for path to change
        page.keyboard.press(f"Key{key}")
        d_1 = line.get_attribute("d")
        expect(line).not_to_have_attribute("d", d_0)

        # Press key again and wait for path to return to original
        page.keyboard.press(f"Key{key}")
        d_2 = line.get_attribute("d")
        expect(line).to_have_attribute("d", d_0)

        # Verify the transformations
        assert d_0 != d_1, "scaling toggle did not change state on first toggle"
        assert d_0 == d_2, "scaling toggle did not return to original state after 2 toggles"

    _test_scaling_toggle_indempotent(spectra_viewer, series=0, key="Q")
    _test_scaling_toggle_indempotent(spectra_viewer, series=1, key="E")


def test_custom_laser_exnorm_persistence(live_server: LiveServer, page: Page) -> None:
    """Test that custom laser exNorm checkbox state persists across page refresh.

    Regression test for bug where "Norm em. to this" checkbox state
    was not saved to sessionStorage.

    Verifies:
    1. Custom laser exNorm checkbox can be checked
    2. Checkbox state persists after page refresh
    """
    create_egfp()

    # Navigate to spectra viewer
    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Add EGFP to the viewer
    _add_egfp_to_viewer(page)

    # Navigate to Light Sources tab
    page.get_by_role("tab", name=re.compile("Light Sources")).click()

    # Click "Add Laser" button
    page.get_by_role("button", name="Add Laser").click()
    page.wait_for_timeout(200)  # Wait for laser to be added

    # Find and check the "Norm em. to this" checkbox
    exnorm_checkbox = page.get_by_role("checkbox", name="Norm em. to this")
    expect(exnorm_checkbox).to_be_visible()
    expect(exnorm_checkbox).not_to_be_checked()

    # Check the checkbox
    exnorm_checkbox.click()
    expect(exnorm_checkbox).to_be_checked()

    # Refresh the page
    page.reload()
    page.wait_for_load_state("networkidle")

    # Verify laser is still there and checkbox is still checked
    page.get_by_role("tab", name=re.compile("Light Sources")).click()
    exnorm_checkbox = page.get_by_role("checkbox", name="Norm em. to this")
    expect(exnorm_checkbox).to_be_visible()
    expect(exnorm_checkbox).to_be_checked()


def test_hidden_spectra_in_share_url(live_server: LiveServer, page: Page) -> None:
    """Test that hidden spectra state is included in shared URLs.

    Regression test for bug where hiddenSpectra was not serialized
    into share URLs.

    This test verifies that hiddenSpectra from the store is properly
    serialized into share URLs when loaded from URL parameters.

    Verifies:
    1. URL with hidden spectra parameter loads correctly
    2. Share URL preserves hidden spectra parameter
    """
    egfp = create_egfp()
    ex_id = str(egfp.default_state.ex_spectrum.id)
    em_id = str(egfp.default_state.em_spectrum.id)

    # Navigate with both spectra, with EX marked as hidden
    url = (
        f"{live_server.url}{reverse('proteins:spectra')}"
        f"?s={ex_id},{em_id}"
        f"&h={ex_id}"  # Mark EX spectrum as hidden via URL parameter
    )
    page.goto(url)
    page.wait_for_load_state("networkidle")

    # Wait for chart to load - should have at least EM spectrum
    chart = page.locator(".highcharts-container")
    expect(chart).to_be_visible()

    # Open share dialog
    page.wait_for_timeout(500)  # Give time for state to settle
    share_button = page.locator('button[aria-controls="simple-menu"]')
    expect(share_button).to_be_visible()
    share_button.click()

    share_menu_item = page.get_by_role("menuitem", name="Share chart as URL")
    expect(share_menu_item).to_be_visible(timeout=10000)
    share_menu_item.click()

    # Get the shared URL
    dialog = page.get_by_role("dialog", name=re.compile("recreate the current graph"))
    expect(dialog).to_be_visible()

    url_textbox = dialog.get_by_role("textbox", name="URL")
    expect(url_textbox).to_be_visible()
    shared_url = url_textbox.input_value()

    # Verify URL contains hidden spectra parameter
    # This is the key test: hiddenSpectra from store should be serialized
    assert "h=" in shared_url, "Shared URL should contain hidden spectra parameter"
    assert ex_id in shared_url, f"Shared URL should contain hidden EX spectrum ID: {ex_id}"


@pytest.mark.parametrize(
    "params",
    [
        {},
        {"xMin": 400, "xMax": 600},
        {"xMin": 400, "xMax": 600, "scaleEC": 1},
    ],
)
def test_spectra_graph(live_server: LiveServer, page: Page, params: dict) -> None:
    egfp = create_egfp()

    ids = f"{egfp.default_state.ex_spectrum.id},{egfp.default_state.em_spectrum.id}"
    query = f"?s={ids}" + "".join(f"&{k}={v}" for k, v in params.items())
    url = f"{live_server.url}{reverse('proteins:spectra_graph')}{query}"
    page.goto(url)
    expect(page).to_have_url(url)

    spectra_viewer = page.locator("#spectra-viewer")
    expect(spectra_viewer).to_be_visible()

    series_paths = spectra_viewer.locator(".highcharts-series-group path.highcharts-area")
    expect(series_paths).to_have_count(2)

    if "xMin" in params and "xMax" in params:
        # Verify that x-axis extremes match the parameters
        x_axis_labels = spectra_viewer.locator(".highcharts-xaxis-labels text")
        label_texts = x_axis_labels.all_text_contents()
        label_values = [int(text) for text in label_texts if text.strip() and text.strip().isdigit()]
        assert min(label_values) >= params["xMin"], f"Min label {min(label_values)} should be >= {params['xMin']}"
        assert max(label_values) <= params["xMax"], f"Max label {max(label_values)} should be <= {params['xMax']}"


def test_subtype_visibility_toggles(spectra_viewer: Page) -> None:
    """Test that clicking subtype toggle buttons (EX, EM, 2P) changes spectrum visibility.

    Regression test for bug where EX and 2P toggle buttons did not respond to clicks
    while EM toggle worked correctly.
    """
    # Add EGFP to the viewer
    _add_egfp_to_viewer(spectra_viewer)

    # Wait for chart to load with spectra
    visible_areas = spectra_viewer.locator(".highcharts-series-group path.highcharts-area:visible")
    expect(visible_areas).to_have_count(2)  # EGFP has EX and EM spectra

    # Find the subtype toggle buttons within the EGFP selector group
    tab_wrapper = spectra_viewer.locator(".tab-wrapper")
    ex_button = tab_wrapper.get_by_role("button", name="EX", exact=True)
    em_button = tab_wrapper.get_by_role("button", name="EM", exact=True)
    twop_button = tab_wrapper.get_by_role("button", name="2P", exact=True)

    # Verify both buttons are visible
    expect(ex_button).to_be_visible()
    expect(em_button).to_be_visible()
    expect(twop_button).to_be_visible()

    # Click EX button to hide EX spectrum
    ex_button.click()
    expect(visible_areas).to_have_count(1)  # Only EM should be visible

    # Click EM button to hide EM spectrum
    em_button.click()
    expect(visible_areas).to_have_count(0)  # No spectra visible

    # Click EX button again to show EX spectrum
    ex_button.click()
    expect(visible_areas).to_have_count(1)  # EX is visible again

    # Click EM button again to show EM spectrum
    em_button.click()
    expect(visible_areas).to_have_count(2)  # Both EX and EM visible again

    # # Click 2P button again to show 2P spectrum
    # twop_button.click()
    # expect(visible_areas).to_have_count(3)  # EX, EM, and 2P visible
    # twop_button.click()
    # expect(visible_areas).to_have_count(2)  # Back to EX and EM


def test_x_range_pickers(spectra_viewer: Page) -> None:
    """Test that X range picker inputs correctly set chart zoom.

    Verifies:
    1. Range picker inputs are visible when spectra are loaded
    2. Typing new min/max values and pressing Enter updates the chart
    3. Chart extremes reflect the entered values
    4. Values persist across page refresh
    5. URL parameters (xMin/xMax) correctly set initial range
    """
    # Add EGFP to the viewer
    _add_egfp_to_viewer(spectra_viewer)

    # Wait for chart to load
    chart = spectra_viewer.locator(".highcharts-container")
    expect(chart).to_be_visible()

    # Range pickers should be visible
    min_input = spectra_viewer.locator('input[name="min"]')
    max_input = spectra_viewer.locator('input[name="max"]')
    expect(min_input).to_be_visible()
    expect(max_input).to_be_visible()

    # Get initial values (should be autoscaled)
    initial_min = min_input.input_value()
    initial_max = max_input.input_value()
    assert initial_min, "Min input should have a value"
    assert initial_max, "Max input should have a value"
    # Set min to 400
    min_input.click()
    min_input.fill("400")
    min_input.press("Enter")
    expect(min_input).to_have_value("400")

    # Set max to 600
    max_input.click()
    max_input.fill("600")
    max_input.press("Enter")
    expect(max_input).to_have_value("600")

    # Verify chart has zoomed by checking x-axis labels
    # Should see labels around 400-600 range, not the original full range
    x_axis_labels = spectra_viewer.locator(".highcharts-xaxis-labels text")
    # Note: Some labels may be hidden by XRangePickers positioning logic
    # Get all label values (including hidden ones) and verify they're in the expected range
    label_texts = x_axis_labels.all_text_contents()
    label_values = [int(text) for text in label_texts if text.strip() and text.strip().isdigit()]
    assert len(label_values) > 0, "Should have numeric x-axis labels"
    assert min(label_values) >= 400, f"Min label {min(label_values)} should be >= 400"
    assert max(label_values) <= 600, f"Max label {max(label_values)} should be <= 600"

    # Test persistence: refresh the page and verify values are restored
    spectra_viewer.reload()
    spectra_viewer.wait_for_load_state("networkidle")

    # Wait for chart to reload
    expect(chart).to_be_visible()
    expect(min_input).to_be_visible()
    expect(max_input).to_be_visible()

    # Verify values persisted
    expect(min_input).to_have_value("400")
    expect(max_input).to_have_value("600")

    # Verify chart is still zoomed after reload
    label_texts_after = x_axis_labels.all_text_contents()
    label_values_after = [int(text) for text in label_texts_after if text.strip() and text.strip().isdigit()]
    min_label_after = min(label_values_after)
    max_label_after = max(label_values_after)
    assert min_label_after >= 400, f"Min label after reload {min_label_after} should be >= 400"
    assert max_label_after <= 600, f"Max label after reload {max_label_after} should be <= 600"

    # click the reset zoom button to return to full range
    zoom_btn = spectra_viewer.locator("g.highcharts-reset-zoom").first
    expect(zoom_btn).to_be_visible()
    zoom_btn.click()
    expect(min_input).to_have_value(initial_min)
    expect(max_input).to_have_value(initial_max)

    # verify chart is back to full range
    label_texts_full = x_axis_labels.all_text_contents()
    label_values_full = [int(text) for text in label_texts_full if text.strip() and text.strip().isdigit()]
    assert min(label_values_full) == 300
    assert max(label_values_full) == 900


def test_x_range_url_parameters(live_server: LiveServer, page: Page) -> None:
    """Test that xMin and xMax URL parameters correctly set initial chart range.

    Verifies:
    1. URL with xMin/xMax parameters loads correctly
    2. Chart range matches URL parameters
    3. Range picker inputs show correct values
    """
    egfp = create_egfp()
    ex_id = egfp.default_state.ex_spectrum.id
    em_id = egfp.default_state.em_spectrum.id

    # Navigate with URL parameters including xMin and xMax
    url = f"{live_server.url}{reverse('proteins:spectra')}?s={ex_id},{em_id}&xMin=420&xMax=580"
    page.goto(url)
    page.wait_for_load_state("networkidle")

    # Wait for chart to load
    chart = page.locator(".highcharts-container")
    expect(chart).to_be_visible()

    # Verify range picker inputs show the URL parameter values
    min_input = page.locator('input[name="min"]')
    max_input = page.locator('input[name="max"]')
    expect(min_input).to_be_visible()
    expect(max_input).to_be_visible()
    expect(min_input).to_have_value("420")
    expect(max_input).to_have_value("580")

    # Verify chart is zoomed by checking x-axis labels
    x_axis_labels = page.locator(".highcharts-xaxis-labels text")

    label_texts = x_axis_labels.all_text_contents()
    label_values = [int(text) for text in label_texts if text.strip() and text.strip().isdigit()]
    assert len(label_values) > 0, "Should have numeric x-axis labels"
    assert min(label_values) >= 420, f"Min label {min(label_values)} should be >= 420"
    assert max(label_values) <= 580, f"Max label {max(label_values)} should be <= 580"


def test_remove_all_spectra(spectra_viewer: Page) -> None:
    """Test that removing all spectra from the viewer works correctly."""
    # Add EGFP to the viewer
    x_axis_labels = spectra_viewer.locator(".highcharts-xaxis-labels")
    zoom_button = spectra_viewer.locator("text").filter(has_text="Reset zoom")
    expect(x_axis_labels).not_to_be_visible()
    expect(zoom_button).not_to_be_visible()

    _add_egfp_to_viewer(spectra_viewer)

    # manually set the xMax
    max_input = spectra_viewer.locator('input[name="max"]')
    max_input.click()
    max_input.fill("600")
    max_input.press("Enter")
    expect(max_input).to_have_value("600")
    expect(x_axis_labels).to_be_visible()
    expect(zoom_button).to_be_visible()

    # blur the input to ensure the value is set
    max_input.blur()

    # open config drawer with ","
    spectra_viewer.keyboard.press("Comma")
    reset_btn = spectra_viewer.locator('button:has-text("Remove All Spectra")')
    expect(reset_btn).to_be_visible()
    reset_btn.click()
    expect(x_axis_labels).not_to_be_visible()
    expect(zoom_button).not_to_be_visible()
    series_paths = spectra_viewer.locator(".highcharts-series-group path.highcharts-area")
    expect(series_paths).to_have_count(0)
