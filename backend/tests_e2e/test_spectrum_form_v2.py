"""End-to-end tests for SpectrumFormV2 (enhanced spectrum submission form)."""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import (
    FilterFactory,
    ProteinFactory,
    SpectrumFactory,
    StateFactory,
)
from proteins.models import Spectrum, State

if TYPE_CHECKING:
    from collections.abc import Iterator

    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


# --- Fixtures ---


@pytest.fixture
def sample_csv_file() -> Iterator[Path]:
    """Create a temporary CSV file with sample spectrum data."""
    content = """wavelength,EGFP_ex,EGFP_em,Filter_BP525
400,0.15,0.00,0.0
450,0.80,0.02,0.0
488,1.00,0.15,0.0
507,0.60,1.00,0.95
530,0.15,0.50,0.85
600,0.00,0.00,0.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        path = Path(f.name)
    yield path
    path.unlink(missing_ok=True)


@pytest.fixture
def european_csv_file() -> Iterator[Path]:
    """European-format CSV (semicolon delimiter, comma decimal)."""
    content = """Wavelength;Relative intensity
476,0599976;11,14449883
494,8699951;1,0
500,2099991;0,143111453
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        path = Path(f.name)
    yield path
    path.unlink(missing_ok=True)


@pytest.fixture
def spectrum_form_page(live_server: LiveServer, auth_page: Page) -> Iterator[Page]:
    """Navigate to the spectrum form page."""
    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page.locator("#spectrum-form-v2")).to_be_visible()
    yield auth_page


# --- Helpers ---


def _upload_csv(page: Page, csv_path: Path) -> None:
    """Upload a CSV file to the form."""
    page.locator("#id_file").set_input_files(str(csv_path))
    expect(page.locator("#column-picker-container")).to_be_visible()


def _select_columns(page: Page, wavelength_col: int, data_cols: list[int]) -> None:
    """Select wavelength and data columns, then continue."""
    headers = page.locator("#column-picker-container th.column-header")
    headers.nth(wavelength_col).click()
    for col in data_cols:
        headers.nth(col).click()
    page.locator("#continue-btn").click()
    expect(page.locator("#spectra-preview-container")).to_be_visible()


def _fill_spectrum_card(
    page: Page,
    category: str,
    subtype: str,
    owner: str,
    card_index: int = 0,
) -> None:
    """Fill out a spectrum card with the given values."""
    card = page.locator(".spectrum-card").nth(card_index)
    card.locator('[id^="category-select-"]').select_option(category)
    card.locator('[id^="subtype-select-"]').select_option(subtype)
    owner_input = card.locator('[id^="owner-input-"]')
    if owner_input.is_visible():
        owner_input.fill(owner)


def _fill_source(page: Page, source: str) -> None:
    """Fill the source field."""
    page.locator("#id_source").fill(source)


# --- Basic Page Tests ---


def test_spectrum_form_page_requires_auth(live_server: LiveServer, page: Page) -> None:
    """Unauthenticated users are redirected to login."""
    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    page.goto(url)
    expect(page).to_have_url(re.compile(r".*/accounts/login/.*"))


def test_spectrum_form_page_loads(spectrum_form_page: Page) -> None:
    """Form page loads with expected elements."""
    page = spectrum_form_page
    expect(page.locator("h2")).to_contain_text("Submit Spectra")
    expect(page.locator("#id_file")).to_be_visible()
    expect(page.locator("#submit-btn")).to_be_disabled()


# --- File Upload and Column Picker ---


def test_file_upload_shows_column_picker(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Uploading a file shows column picker with correct headers."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)

    column_picker = page.locator("#column-picker-container")
    expect(column_picker).to_contain_text("Select columns")
    expect(column_picker.locator("th.column-header")).to_have_count(4)
    expect(page.locator("#continue-btn")).to_be_disabled()


def test_column_selection_workflow(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Column selection enables continue button and shows spectrum cards."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    headers = page.locator("#column-picker-container th.column-header")

    # Select wavelength column
    headers.nth(0).click()
    expect(page.locator("#wave-status")).to_contain_text("wavelength")
    expect(page.locator("#continue-btn")).to_be_disabled()

    # Select data column enables continue
    headers.nth(1).click()
    expect(page.locator("#continue-btn")).to_be_enabled()

    # Continue shows spectrum cards
    page.locator("#continue-btn").click()
    expect(page.locator(".spectrum-card")).to_have_count(1)


def test_multiple_column_selection(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Selecting multiple data columns creates multiple spectrum cards."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1, 2, 3])
    expect(page.locator(".spectrum-card")).to_have_count(3)


def test_european_csv_format_parsing(spectrum_form_page: Page, european_csv_file: Path) -> None:
    """European CSV format (semicolon delimiter, comma decimal) parses correctly."""
    page = spectrum_form_page
    _upload_csv(page, european_csv_file)

    headers = page.locator("#column-picker-container th.column-header")
    expect(headers).to_have_count(2)

    # Verify numeric values are parsed (not "-" which indicates NaN)
    first_cell = page.locator("#column-picker-container table tbody tr").first.locator("td").first
    assert "476." in (first_cell.text_content() or "")

    _select_columns(page, wavelength_col=0, data_cols=[1])
    expect(page.locator(".spectrum-card")).to_have_count(1)


# --- Spectrum Card Category/Subtype Tests ---


@pytest.mark.parametrize(
    ("category", "expected_subtype", "subtype_count"),
    [
        ("c", "qe", 1),  # Camera auto-selects only subtype
        ("l", "pd", 1),  # Light auto-selects only subtype
    ],
)
def test_category_auto_selects_single_subtype(
    spectrum_form_page: Page,
    sample_csv_file: Path,
    category: str,
    expected_subtype: str,
    subtype_count: int,
) -> None:
    """Categories with single subtype auto-select it."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option(category)
    subtype_select = card.locator('[id^="subtype-select-"]')

    expect(subtype_select).to_have_value(expected_subtype)
    expect(subtype_select.locator("option")).to_have_count(subtype_count)


@pytest.mark.parametrize(
    ("category", "bio_fields_visible"),
    [
        ("p", True),  # Protein - bio fields visible
        ("d", True),  # Dye - bio fields visible
        ("f", False),  # Filter - no bio fields
        ("c", False),  # Camera - no bio fields
    ],
)
def test_bio_fields_visibility(
    spectrum_form_page: Page,
    sample_csv_file: Path,
    category: str,
    bio_fields_visible: bool,
) -> None:
    """pH and Solvent fields show only for bio categories (Protein, Dye)."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option(category)

    ph_input = card.locator('[id^="ph-input-"]')
    if bio_fields_visible:
        expect(ph_input).to_be_visible()
    else:
        expect(ph_input).not_to_be_visible()


@pytest.mark.parametrize(
    ("category", "peak_visible"),
    [
        ("p", True),  # Protein - normalized
        ("l", True),  # Light - normalized
        ("f", False),  # Filter - not normalized
        ("c", False),  # Camera - not normalized
    ],
)
def test_peak_badge_visibility(
    spectrum_form_page: Page,
    sample_csv_file: Path,
    category: str,
    peak_visible: bool,
) -> None:
    """Peak badge shows only for normalizable categories."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option(category)
    peak_badge = card.locator('[id^="peak-badge-"]')

    if peak_visible:
        expect(peak_badge).to_be_visible()
    else:
        expect(peak_badge).not_to_be_visible()


def test_owner_field_autocomplete_for_protein(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Protein category uses Select2 autocomplete for owner."""
    page = spectrum_form_page
    ProteinFactory(name="TestProtein", slug="testprotein")

    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("p")

    expect(card.locator('[id^="owner-input-"]')).not_to_be_visible()
    expect(card.locator(".select2-container")).to_be_visible()


def test_protein_category_form_submits_without_focus_error(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Protein category with hidden owner-input doesn't cause 'not focusable' error on submit.

    Regression test: When protein category is selected, the text input is hidden and
    a Select2 dropdown is shown. The hidden input must not have 'required' attribute,
    otherwise the browser fails to submit with 'An invalid form control is not focusable'.
    """
    page = spectrum_form_page
    protein = ProteinFactory(name="SubmitTestProtein", slug="submittestprotein")

    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("p")
    card.locator('[id^="subtype-select-"]').select_option("ex")

    # Select protein using Select2
    card.locator(".select2-container").click()
    page.locator(".select2-search__field").fill("SubmitTest")
    page.locator(".select2-results__option").first.click()

    # Verify the hidden input does NOT have required attribute
    owner_input = card.locator('[id^="owner-input-"]')
    expect(owner_input).not_to_be_visible()
    expect(owner_input).not_to_have_attribute("required", "")

    # Fill source and submit
    _fill_source(page, "Protein submit test")
    page.locator("#id_confirmation").check()
    page.locator("#submit-btn").click()

    # Should redirect successfully (no focus error)
    expect(page).to_have_url(re.compile(r".*/spectra/submitted/"))

    # Verify spectrum was created
    spectrum = Spectrum.objects.all_objects().filter(source="Protein submit test").first()
    assert spectrum is not None
    assert spectrum.category == "p"
    assert spectrum.owner_fluor is not None
    # owner_fluor is a State (subclass of FluorState) for proteins
    state = State.objects.get(pk=spectrum.owner_fluor.pk)
    assert state.protein.slug == protein.slug


# --- Validation Tests ---


def test_validation_messages(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Validation shows appropriate messages as fields are filled."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    validation_msg = page.locator("#validation-message")
    card = page.locator(".spectrum-card").first

    # Missing category
    expect(validation_msg).to_contain_text("category")

    # Fill category, now missing subtype/owner
    card.locator('[id^="category-select-"]').select_option("f")
    card.locator('[id^="subtype-select-"]').select_option("bp")

    # Missing owner
    expect(validation_msg).to_contain_text("owner")

    # Fill owner, now missing source
    card.locator('[id^="owner-input-"]').fill("Test Filter")
    expect(validation_msg).to_contain_text("Source")

    # Fill source, now missing confirmation
    page.locator("#id_source").fill("Test source")
    expect(validation_msg).to_contain_text("confirmation")

    # Check confirmation - ready to submit
    page.locator("#id_confirmation").check()
    expect(validation_msg).to_contain_text("ready to submit")


def test_doi_validation(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """DOI validation on reference field."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])
    _fill_spectrum_card(page, category="f", subtype="bp", owner="Test Filter")

    reference_input = page.locator("#id_reference")
    validation_msg = page.locator("#validation-message")

    # Invalid DOI
    reference_input.fill("not-a-doi")
    reference_input.blur()
    expect(validation_msg).to_contain_text("DOI")

    # Valid DOI - still needs confirmation
    reference_input.fill("10.1234/test.doi")
    reference_input.blur()
    expect(validation_msg).to_contain_text("confirmation")

    # Check confirmation - ready to submit
    page.locator("#id_confirmation").check()
    expect(validation_msg).to_contain_text("ready to submit")


def test_submit_button_enables_when_valid(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Submit button enables when all required fields are valid."""
    page = spectrum_form_page
    submit_btn = page.locator("#submit-btn")

    expect(submit_btn).to_be_disabled()

    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])
    expect(submit_btn).to_be_disabled()

    _fill_spectrum_card(page, category="f", subtype="bp", owner="Test Filter")
    expect(submit_btn).to_be_disabled()

    _fill_source(page, "Test source")
    expect(submit_btn).to_be_disabled()  # Still needs confirmation

    page.locator("#id_confirmation").check()
    expect(submit_btn).to_be_enabled()


# --- Card Actions ---


def test_remove_spectrum(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Removing a spectrum card updates the count."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1, 2])

    expect(page.locator(".spectrum-card")).to_have_count(2)

    page.locator(".spectrum-card").first.locator('[id^="remove-btn-"]').click()
    expect(page.locator(".spectrum-card")).to_have_count(1)


def test_remove_all_spectra_resets_to_column_picker(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Removing all spectra resets UI to column picker and allows re-upload."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    # Verify we're in spectrum card view
    expect(page.locator("#global-source-fields")).to_be_visible()
    expect(page.locator("#column-picker-container")).not_to_be_visible()

    # Remove the only card
    page.locator('[id^="remove-btn-"]').first.click()

    # Should reset to column picker view
    expect(page.locator("#global-source-fields")).not_to_be_visible()
    expect(page.locator("#spectra-preview-container")).not_to_be_visible()
    expect(page.locator("#column-picker-container")).to_be_visible()

    # Should be able to re-upload the same file
    _upload_csv(page, sample_csv_file)
    expect(page.locator("#column-picker-container th.column-header")).to_have_count(4)


def test_status_indicators_update(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Status indicators update as fields are filled."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    status_icon = card.locator('[id^="status-icon-"]')

    expect(status_icon).to_contain_text("⚠️")

    _fill_spectrum_card(page, category="f", subtype="bp", owner="Test Filter")
    expect(status_icon).to_contain_text("✅")


def test_chart_click_sets_peak_wavelength(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Clicking chart updates peak wavelength.

    For normalizable categories (dye, protein, light), clicking on the chart
    sets a manual peak wavelength. The peak badge should update to reflect
    the new peak location.
    """
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])  # EGFP_ex column

    card = page.locator(".spectrum-card").first
    peak_badge = card.locator('[id^="peak-badge-"]')

    # Select dye category (normalizable, shows peak badge)
    card.locator('[id^="category-select-"]').select_option("d")
    card.locator('[id^="subtype-select-"]').select_option("ex")

    # Initial peak should be at 488nm (max value in sample data)
    expect(peak_badge).to_be_visible()
    expect(peak_badge).to_contain_text("Peak: 488 nm")

    # Click on the left side of the chart to set a new peak
    # The chart spans wavelengths 400-600nm
    chart_container = card.locator('[id^="chart-container-"]')
    box = chart_container.bounding_box()
    assert box is not None

    # Click at ~20% from left edge (should be around 440nm area)
    click_x = box["x"] + box["width"] * 0.2
    click_y = box["y"] + box["height"] / 2
    page.mouse.click(click_x, click_y)

    # Peak should change from the initial 488nm to something in the 400-480nm range
    # The exact value depends on chart rendering, so we just verify it changed
    expect(peak_badge).not_to_contain_text("Peak: 488 nm")
    # Verify it's showing a valid peak (not "No peak" or "--")
    expect(peak_badge).to_contain_text(re.compile(r"Peak: \d+ nm"))


# --- Form Submission Tests ---


def test_submit_filter_spectrum_with_source(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Submit a filter spectrum with source field (no reference)."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[3])  # Filter column

    _fill_spectrum_card(page, category="f", subtype="bp", owner="E2E Test Filter")
    _fill_source(page, "E2E Test Source")

    # Check confirmation and submit
    page.locator("#id_confirmation").check()
    page.locator("#submit-btn").click()

    # Should redirect to success page
    expect(page).to_have_url(re.compile(r".*/spectra/submitted/"))
    expect(page.locator("body")).to_contain_text("Thank You!")

    # Verify spectrum was created in database
    spectrum = Spectrum.objects.all_objects().filter(source="E2E Test Source").first()
    assert spectrum is not None
    assert spectrum.category == "f"
    assert spectrum.subtype == "bp"
    assert spectrum.reference is None  # No reference was provided


def test_submit_light_spectrum(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Submit a light source spectrum (exercises light category handling)."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[3])

    _fill_spectrum_card(page, category="l", subtype="pd", owner="E2E Test Light")
    _fill_source(page, "E2E Test Light Source")

    # Check confirmation and submit
    page.locator("#id_confirmation").check()
    expect(page.locator("#submit-btn")).to_be_enabled()
    page.locator("#submit-btn").click()

    # Should redirect to success page
    expect(page).to_have_url(re.compile(r".*/spectra/submitted/"))

    # Verify spectrum was created
    spectrum = Spectrum.objects.all_objects().filter(owner_light__name="E2E Test Light").first()
    assert spectrum is not None
    assert spectrum.category == "l"
    assert spectrum.subtype == "pd"


def test_submit_multiple_spectra(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Submit multiple spectra at once."""
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1, 3])  # 2 columns

    # Fill first card (will be a protein excitation - needs protein owner)
    # Use filter instead for simpler test
    _fill_spectrum_card(page, category="f", subtype="bp", owner="E2E Multi Filter 1", card_index=0)
    _fill_spectrum_card(page, category="f", subtype="lp", owner="E2E Multi Filter 2", card_index=1)
    _fill_source(page, "E2E Multi Test")

    page.locator("#id_confirmation").check()
    page.locator("#submit-btn").click()

    expect(page).to_have_url(re.compile(r".*/spectra/submitted/"))

    # Verify both spectra were created
    spectra = Spectrum.objects.all_objects().filter(source="E2E Multi Test")
    assert spectra.count() == 2


# --- Peak Selection Tests ---


@pytest.fixture
def peak_snap_csv_file() -> Iterator[Path]:
    """CSV with clear, sharp peak for testing peak snap behavior.

    Data has a clear peak at 500nm. When clicking near but not on the peak,
    the displayed peak should "snap" to the local maximum (500nm).
    """
    # Create a gaussian-like peak centered at 500nm
    content = "wavelength,intensity\n"
    for wave in range(400, 601, 5):
        # Gaussian centered at 500nm with sigma=30
        value = 0.1 + 0.9 * (2.718 ** (-(((wave - 500) / 30) ** 2)))
        content += f"{wave},{value:.4f}\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        path = Path(f.name)
    yield path
    path.unlink(missing_ok=True)


def test_peak_snaps_to_local_maximum(spectrum_form_page: Page, peak_snap_csv_file: Path) -> None:
    """Clicking near a peak should snap the displayed peak to the local maximum.

    This is critical UX behavior: when a user clicks on the chart to set a peak,
    we find the local maximum within ±25nm of the click position. The displayed
    peak wavelength and the normalization target must both use this local maximum,
    not the raw click coordinates.

    The normalized value at the displayed peak should always be exactly 1.0.
    """
    page = spectrum_form_page
    _upload_csv(page, peak_snap_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    peak_badge = card.locator('[id^="peak-badge-"]')

    # Select dye category (normalizable) and fill required fields
    card.locator('[id^="category-select-"]').select_option("d")
    card.locator('[id^="subtype-select-"]').select_option("ex")
    card.locator('[id^="owner-input-"]').fill("Test Dye")

    # Initial peak should be around 500nm (may vary slightly due to interpolation)
    expect(peak_badge).to_be_visible()
    initial_badge_text = peak_badge.text_content()
    assert initial_badge_text is not None
    initial_peak_match = re.search(r"Peak: (\d+) nm", initial_badge_text)
    assert initial_peak_match, f"Could not parse peak from badge: {initial_badge_text}"
    initial_peak = int(initial_peak_match.group(1))
    assert 498 <= initial_peak <= 502, f"Initial peak should be ~500nm, got {initial_peak}nm"

    # Click at right side of chart (high wavelength region, ~560-580nm)
    # The chart spans 400-600nm, so clicking at 80% from left = ~560nm
    chart_container = card.locator('[id^="chart-container-"]')
    box = chart_container.bounding_box()
    assert box is not None

    # Click at 80% from left edge - should be around 560nm
    # The local max within ±25nm of 560nm is still near 500nm (gaussian peak)
    click_x = box["x"] + box["width"] * 0.8
    click_y = box["y"] + box["height"] / 2
    page.mouse.click(click_x, click_y)

    # Get the new peak after click
    new_badge_text = peak_badge.text_content()
    assert new_badge_text is not None
    new_peak_match = re.search(r"Peak: (\d+) nm", new_badge_text)
    assert new_peak_match, f"Could not parse peak from badge: {new_badge_text}"
    new_peak = int(new_peak_match.group(1))

    # CRITICAL: If clicking at ~560nm, the raw click would be way off from 500nm.
    # But the local max search (±25nm) should find a peak in 535-585nm range.
    # For a gaussian centered at 500nm, the local max in that range would be ~535nm.
    # The key point: it should NOT be the raw click position (560nm), but rather
    # a local maximum that represents an actual peak in the data.

    # Verify the JSON data has the correct peak and normalized value
    spectra_json = page.locator("#id_spectra_json").input_value()
    assert spectra_json, "spectra_json should not be empty"

    spectra_data = json.loads(spectra_json)
    assert len(spectra_data) == 1

    spectrum = spectra_data[0]
    json_peak = spectrum["peak_wave"]

    # The displayed peak and JSON peak should match
    assert json_peak == new_peak, f"Badge shows {new_peak}nm but JSON has {json_peak}nm"

    # MOST CRITICAL: The normalized value at the displayed peak must be exactly 1.0
    # This proves the displayed peak is actually what was used for normalization
    data = spectrum["data"]
    peak_values = [val for wave, val in data if wave == json_peak]
    assert len(peak_values) == 1, f"No data point at peak wavelength {json_peak}nm"
    assert peak_values[0] == 1.0, (
        f"Expected normalized peak value of 1.0 at {json_peak}nm, got {peak_values[0]}. "
        "This means the displayed peak doesn't match the normalization target!"
    )


def test_peak_snap_to_different_local_max(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Clicking in a different region finds that region's local maximum.

    The sample data has:
    - 450nm: 0.80 (local max in lower wavelength region)
    - 488nm: 1.00 (global max)

    Clicking near 450nm should snap to a local max in that region, not 488nm.
    The exact wavelength may vary due to interpolation, but the key is that
    the normalized value at the displayed peak should be exactly 1.0.
    """
    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])  # EGFP_ex

    card = page.locator(".spectrum-card").first
    peak_badge = card.locator('[id^="peak-badge-"]')

    card.locator('[id^="category-select-"]').select_option("d")
    card.locator('[id^="subtype-select-"]').select_option("ex")
    card.locator('[id^="owner-input-"]').fill("Test Dye")

    # Initial peak should be at 488nm (global max)
    expect(peak_badge).to_contain_text("Peak: 488 nm")

    # Click in the 440-460nm region (around 450nm local max)
    chart_container = card.locator('[id^="chart-container-"]')
    box = chart_container.bounding_box()
    assert box is not None

    # Chart spans 400-600nm (200nm range), click at ~25% from left edge
    # which should be around 450nm area
    click_x = box["x"] + box["width"] * 0.25
    click_y = box["y"] + box["height"] / 2
    page.mouse.click(click_x, click_y)

    # Get the new peak value
    new_badge_text = peak_badge.text_content()
    assert new_badge_text is not None
    new_peak_match = re.search(r"Peak: (\d+) nm", new_badge_text)
    assert new_peak_match, f"Could not parse peak from badge: {new_badge_text}"
    new_peak = int(new_peak_match.group(1))

    # The new peak should be in the 425-475nm range (±25nm from click at ~450nm)
    assert 425 <= new_peak <= 475, (
        f"Expected peak in 425-475nm range, got {new_peak}nm. "
        "Peak should have snapped to local max near click, not global max at 488nm."
    )

    # Verify the JSON data has the correct normalized value at peak
    spectra_json = page.locator("#id_spectra_json").input_value()
    assert spectra_json, "spectra_json should not be empty"

    spectra_data = json.loads(spectra_json)
    spectrum = spectra_data[0]
    json_peak = spectrum["peak_wave"]

    # The displayed peak and JSON peak should match
    assert json_peak == new_peak, f"Badge shows {new_peak}nm but JSON has {json_peak}nm"

    # CRITICAL: The normalized value at the displayed peak must be exactly 1.0
    data = spectrum["data"]
    peak_values = [val for wave, val in data if wave == json_peak]
    assert len(peak_values) == 1, f"No data point at peak wavelength {json_peak}nm"
    assert peak_values[0] == 1.0, (
        f"Expected normalized peak value of 1.0 at {json_peak}nm, got {peak_values[0]}. "
        "This means the displayed peak doesn't match the normalization target!"
    )


# --- Duplicate Detection Tests ---


def test_protein_existing_subtype_is_disabled(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """For proteins, existing subtypes should be disabled in the dropdown.

    When a protein already has an excitation spectrum in the database,
    the "Excitation" option in the subtype dropdown should be disabled
    to prevent submitting duplicate spectra.
    """
    # Create a protein with an existing excitation spectrum
    protein = ProteinFactory(
        name="DuplicateTestProtein", slug="duplicatetestprotein", default_state=None
    )
    state = StateFactory(protein=protein, name="default")
    # Clean up any existing spectra from previous test runs (--reuse-db)
    state.spectra.all().delete()
    SpectrumFactory(category="p", subtype="ex", owner_fluor=state)

    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("p")

    # Select protein using Select2 BEFORE selecting subtype
    card.locator(".select2-container").click()
    page.locator(".select2-search__field").fill("DuplicateTest")
    page.locator(".select2-results__option").first.click()

    # The excitation subtype should now be disabled
    subtype_select = card.locator('[id^="subtype-select-"]')
    ex_option = subtype_select.locator('option[value="ex"]')
    expect(ex_option).to_be_disabled()
    expect(ex_option).to_contain_text("(exists)")

    # Other subtypes should still be enabled
    em_option = subtype_select.locator('option[value="em"]')
    expect(em_option).not_to_be_disabled()

    # For proteins, no warning should show on exact match (user selected from autocomplete,
    # it's expected to add new subtypes to existing proteins)
    owner_warning = card.locator('[id^="owner-warning-"]')
    expect(owner_warning).to_be_hidden()

    # Select an available subtype (emission)
    subtype_select.select_option("em")

    # Fill source and confirmation
    _fill_source(page, "Duplicate test source")
    page.locator("#id_confirmation").check()

    # Submit button should be enabled (we selected an available subtype)
    expect(page.locator("#submit-btn")).to_be_enabled()


def test_protein_clears_invalid_subtype_on_owner_change(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """If user selected a subtype before owner, and that subtype exists, it gets cleared.

    When a user selects excitation subtype first, then selects a protein that already
    has an excitation spectrum, the subtype should be cleared (since it's now disabled).
    """
    # Create a protein with an existing excitation spectrum
    protein = ProteinFactory(
        name="SubtypeTestProtein", slug="subtypetestprotein", default_state=None
    )
    state = StateFactory(protein=protein, name="default")
    # Clean up any existing spectra from previous test runs (--reuse-db)
    state.spectra.all().delete()
    SpectrumFactory(category="p", subtype="ex", owner_fluor=state)

    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("p")

    # Select subtype BEFORE selecting protein
    subtype_select = card.locator('[id^="subtype-select-"]')
    subtype_select.select_option("ex")
    expect(subtype_select).to_have_value("ex")

    # Now select protein using Select2
    card.locator(".select2-container").click()
    page.locator(".select2-search__field").fill("SubtypeTest")
    page.locator(".select2-results__option").first.click()

    # The subtype should be cleared because "ex" is now disabled
    expect(subtype_select).to_have_value("")

    # The excitation option should be disabled
    ex_option = subtype_select.locator('option[value="ex"]')
    expect(ex_option).to_be_disabled()


def test_dye_existing_subtype_is_disabled(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """For dyes, existing subtypes should be disabled in the dropdown.

    When a dye already has an emission spectrum in the database,
    the "Emission" option in the subtype dropdown should be disabled.
    """
    from proteins.models.dye import Dye, DyeState
    from references.factories import ReferenceFactory

    # Use a unique name and clean up any conflicting data from previous runs
    dye_name = "DyeSubtypeTestUnique"
    # Delete ALL DyeStates with this owner_name (not just via Dye cascade)
    # to handle orphaned states from previous test runs
    DyeState.objects.filter(owner_name=dye_name).delete()
    Dye.objects.filter(name=dye_name).delete()

    # Create dye and state directly (not via factory which auto-creates spectra)
    dye = Dye.objects.create(
        name=dye_name,
        slug=dye_name.lower(),
        primary_reference=ReferenceFactory(),
    )
    # Create DyeState directly with name='default' so label equals dye name
    # Don't use DyeStateFactory as it auto-creates ex/em/2p spectra via RelatedFactory
    dye_state = DyeState.objects.create(
        dye=dye,
        name="default",
        slug=f"{dye.slug}_default",
        owner_name=dye_name,
        owner_slug=dye.slug,
    )
    dye.default_state = dye_state
    dye.save()

    # Create only the emission spectrum (don't use SpectrumFactory as it requires em_max)
    from proteins.models import Spectrum

    Spectrum.objects.create(
        category="d",
        subtype="em",
        owner_fluor=dye_state,
        data=[[400, 0.1], [500, 1.0], [600, 0.1]],
    )

    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("d")

    # Enter the exact dye name
    owner_input = card.locator('[id^="owner-input-"]')
    owner_input.fill(dye_name)
    owner_input.blur()  # Trigger validation

    # Wait for the warning to appear (indicates AJAX completed)
    owner_warning = card.locator('[id^="owner-warning-"]')
    expect(owner_warning).to_be_visible()
    expect(owner_warning).to_have_class(re.compile(r"alert-warning"))

    # The emission subtype should now be disabled
    subtype_select = card.locator('[id^="subtype-select-"]')
    em_option = subtype_select.locator('option[value="em"]')
    expect(em_option).to_be_disabled()
    expect(em_option).to_contain_text("(exists)")

    # Other subtypes should still be enabled
    ex_option = subtype_select.locator('option[value="ex"]')
    expect(ex_option).not_to_be_disabled()

    # Select an available subtype (excitation)
    subtype_select.select_option("ex")

    # Fill source and confirmation
    _fill_source(page, "Dye subtype test source")
    page.locator("#id_confirmation").check()

    # Submit button should be enabled
    expect(page.locator("#submit-btn")).to_be_enabled()


def test_filter_exact_match_blocks_submission(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """For filters, exact name match shows error and blocks submission.

    Filters can only have one spectrum per owner, so any exact name match
    should block submission (unlike fluorophores which can have multiple subtypes).
    """

    # Create a filter with an existing spectrum
    FilterFactory(name="ExistingFilter")

    page = spectrum_form_page
    _upload_csv(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    card.locator('[id^="category-select-"]').select_option("f")
    card.locator('[id^="subtype-select-"]').select_option("bp")

    # Enter the exact filter name
    owner_input = card.locator('[id^="owner-input-"]')
    owner_input.fill("ExistingFilter")
    owner_input.blur()  # Trigger validation

    # Should show error (red alert) - exact match found
    owner_warning = card.locator('[id^="owner-warning-"]')
    expect(owner_warning).to_be_visible()
    expect(owner_warning).to_have_class(re.compile(r"alert-danger"))
    expect(owner_warning).to_contain_text("Exact match found")

    # Fill source and confirmation
    _fill_source(page, "Filter test source")
    page.locator("#id_confirmation").check()

    # Submit button should remain disabled due to exact match
    expect(page.locator("#submit-btn")).to_be_disabled()
