"""End-to-end tests for SpectrumFormV2 (enhanced spectrum submission form)."""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import ProteinFactory

if TYPE_CHECKING:
    from collections.abc import Iterator

    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


@pytest.fixture
def sample_csv_file() -> Iterator[Path]:
    """Create a temporary CSV file with sample spectrum data."""
    content = """wavelength,EGFP_ex,EGFP_em,Filter_BP525
400,0.15,0.00,0.0
410,0.25,0.00,0.0
420,0.35,0.00,0.0
430,0.50,0.00,0.0
440,0.65,0.01,0.0
450,0.80,0.02,0.0
460,0.90,0.03,0.0
470,0.95,0.05,0.0
480,0.98,0.08,0.0
488,1.00,0.15,0.0
500,0.85,0.45,0.9
507,0.60,1.00,0.95
510,0.50,0.95,0.92
520,0.30,0.75,0.90
530,0.15,0.50,0.85
540,0.08,0.30,0.0
550,0.04,0.15,0.0
560,0.02,0.08,0.0
570,0.01,0.04,0.0
580,0.00,0.02,0.0
590,0.00,0.01,0.0
600,0.00,0.00,0.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        path = Path(f.name)
    yield path
    path.unlink(missing_ok=True)


@pytest.fixture
def spectrum_form_page(
    live_server: LiveServer, auth_page: Page, sample_csv_file: Path
) -> Iterator[Page]:
    """Navigate to the spectrum form page and upload a sample file."""
    url = f"{live_server.url}{reverse('proteins:submit-spectra-v2')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Wait for the form to be visible
    form = auth_page.locator("#spectrum-form-v2")
    expect(form).to_be_visible()

    yield auth_page
    auth_page.wait_for_load_state("networkidle")


def _upload_csv_file(page: Page, csv_path: Path) -> None:
    """Upload a CSV file to the form."""
    file_input = page.locator("#id_file")
    file_input.set_input_files(str(csv_path))

    # Wait for column picker to appear
    column_picker = page.locator("#column-picker-container")
    expect(column_picker).to_be_visible()


def _select_columns(page: Page, wavelength_col: int, data_cols: list[int]) -> None:
    """Select wavelength column and data columns in the column picker."""
    # Click wavelength column header
    headers = page.locator("#column-picker-container th.column-header")
    headers.nth(wavelength_col).click()

    # Click each data column header
    for col in data_cols:
        headers.nth(col).click()

    # Click continue button
    continue_btn = page.locator("#continue-btn")
    expect(continue_btn).to_be_enabled()
    continue_btn.click()

    # Wait for spectrum cards to appear
    expect(page.locator("#spectra-preview-container")).to_be_visible()


# --- Basic Page Tests ---


def test_spectrum_form_page_requires_auth(live_server: LiveServer, page: Page) -> None:
    """Test that the spectrum form page requires authentication."""
    url = f"{live_server.url}{reverse('proteins:submit-spectra-v2')}"
    page.goto(url)

    # Should redirect to login page
    expect(page).to_have_url(re.compile(r".*/accounts/login/.*"))


def test_spectrum_form_page_loads(spectrum_form_page: Page) -> None:
    """Test that the spectrum form page loads correctly for authenticated users."""
    page = spectrum_form_page

    # Check main elements are present
    expect(page.locator("h2")).to_contain_text("Submit Spectra")
    expect(page.locator("#id_file")).to_be_visible()
    expect(page.locator("#submit-btn")).to_be_visible()
    expect(page.locator("#submit-btn")).to_be_disabled()


# --- File Upload and Column Picker Tests ---


def test_file_upload_shows_column_picker(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that uploading a file shows the column picker."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)

    # Column picker should be visible with instructions
    column_picker = page.locator("#column-picker-container")
    expect(column_picker).to_contain_text("Select columns")

    # Should show all column headers
    headers = column_picker.locator("th.column-header")
    expect(headers).to_have_count(4)  # wavelength + 3 data columns

    # Continue button should be disabled until columns are selected
    expect(page.locator("#continue-btn")).to_be_disabled()


def test_column_selection_workflow(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test the column selection workflow."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)

    headers = page.locator("#column-picker-container th.column-header")

    # Click wavelength column (first column)
    headers.nth(0).click()

    # Status should update
    wave_status = page.locator("#wave-status")
    expect(wave_status).to_contain_text("wavelength")

    # Continue still disabled (no data column selected)
    expect(page.locator("#continue-btn")).to_be_disabled()

    # Click a data column
    headers.nth(1).click()

    # Continue should now be enabled
    expect(page.locator("#continue-btn")).to_be_enabled()

    # Click continue
    page.locator("#continue-btn").click()

    # Spectrum cards should appear
    expect(page.locator(".spectrum-card")).to_have_count(1)


def test_multiple_column_selection(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test selecting multiple data columns."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1, 2, 3])

    # Should have 3 spectrum cards
    expect(page.locator(".spectrum-card")).to_have_count(3)


# --- Spectrum Card Tests ---


def test_spectrum_card_category_selection(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test category selection on spectrum cards."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')

    # Initially should show placeholder
    expect(category_select).to_have_value("")
    expect(subtype_select).to_contain_text("Select category first")

    # Select Protein category
    category_select.select_option("p")
    expect(category_select).to_have_value("p")

    # Subtype should now have options with placeholder
    expect(subtype_select).to_have_value("")
    expect(subtype_select).to_contain_text("-----")

    # Should have protein subtypes available
    expect(subtype_select.locator("option")).to_have_count(5)  # placeholder + 4 subtypes


def test_category_auto_selects_single_subtype(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Test that categories with single subtype auto-select it."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')

    # Select Camera category (has only "qe" subtype)
    category_select.select_option("c")

    # Subtype should be auto-selected
    expect(subtype_select).to_have_value("qe")
    expect(subtype_select.locator("option")).to_have_count(1)

    # Similarly for Light Source
    category_select.select_option("l")
    expect(subtype_select).to_have_value("pd")
    expect(subtype_select.locator("option")).to_have_count(1)


def test_bio_fields_visibility(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that pH and Solvent fields show only for bio categories (Protein, Dye)."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    ph_input = card.locator('[id^="ph-input-"]')
    solvent_input = card.locator('[id^="solvent-input-"]')

    # Select Protein - bio fields should be visible
    category_select.select_option("p")
    expect(ph_input).to_be_visible()
    expect(solvent_input).to_be_visible()

    # Select Filter - bio fields should be hidden
    category_select.select_option("f")
    expect(ph_input).not_to_be_visible()
    expect(solvent_input).not_to_be_visible()

    # Select Dye - bio fields should be visible again
    category_select.select_option("d")
    expect(ph_input).to_be_visible()
    expect(solvent_input).to_be_visible()


def test_peak_badge_visibility(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that peak badge shows only for normalizable categories."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    peak_badge = card.locator('[id^="peak-badge-"]')

    # Select Protein - peak badge should be visible
    category_select.select_option("p")
    expect(peak_badge).to_be_visible()

    # Select Filter - peak badge should be hidden (no normalization)
    category_select.select_option("f")
    expect(peak_badge).not_to_be_visible()

    # Select Camera - peak badge should be hidden (no normalization)
    category_select.select_option("c")
    expect(peak_badge).not_to_be_visible()

    # Select Light - peak badge should be visible (lights are normalized)
    category_select.select_option("l")
    expect(peak_badge).to_be_visible()


def test_scale_factor_units(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that scale factor units change based on subtype."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    units_label = card.locator('[id^="scale-factor-units-"]')

    # Select Protein with Excitation subtype
    category_select.select_option("p")
    subtype_select.select_option("ex")
    expect(units_label).to_contain_text("EC (M⁻¹cm⁻¹)")

    # Change to Emission subtype
    subtype_select.select_option("em")
    expect(units_label).to_contain_text("QE (0-1)")

    # Change to Two-Photon subtype
    subtype_select.select_option("2p")
    expect(units_label).to_contain_text("Cross Section (GM)")


def test_owner_field_text_input_for_non_protein(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Test that non-protein categories use text input for owner."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')
    owner_select = card.locator('[id^="owner-select-"]')

    # Select Filter category
    category_select.select_option("f")

    # Text input should be visible, select should be hidden
    expect(owner_input).to_be_visible()
    expect(owner_select).not_to_be_visible()

    # Can type in the input
    owner_input.fill("Test Filter BP525")
    expect(owner_input).to_have_value("Test Filter BP525")


def test_owner_field_autocomplete_for_protein(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Test that Protein category uses Select2 autocomplete for owner."""
    page = spectrum_form_page

    # Create a protein for autocomplete to find
    ProteinFactory(name="TestProtein", slug="testprotein")

    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')
    card.locator('[id^="owner-select-"]')

    # Select Protein category
    category_select.select_option("p")

    # Text input should be hidden, select should be visible
    expect(owner_input).not_to_be_visible()

    # Select2 should initialize - look for the container
    select2_container = card.locator(".select2-container")
    expect(select2_container).to_be_visible()


# --- Validation Tests ---


def test_validation_requires_category(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that validation requires category selection."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    # Check validation message mentions category
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("category")


def test_validation_requires_owner(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that validation requires owner."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')

    # Select category and subtype but not owner
    category_select.select_option("f")
    subtype_select.select_option("bp")

    # Validation should mention owner
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("owner")


def test_validation_requires_source_or_reference(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Test that validation requires at least one of source or reference."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')

    # Fill all spectrum fields
    category_select.select_option("f")
    subtype_select.select_option("bp")
    owner_input.fill("Test Filter")

    # Validation should mention source/reference
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("Source")
    expect(validation_msg).to_contain_text("Primary Reference")

    # Fill source
    source_input = page.locator("#id_source")
    source_input.fill("Test source")

    # Should now show ready message
    expect(validation_msg).to_contain_text("ready to submit")


def test_doi_validation(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test DOI validation on reference field."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')

    # Fill spectrum fields
    category_select.select_option("f")
    subtype_select.select_option("bp")
    owner_input.fill("Test Filter")

    # Enter invalid DOI
    reference_input = page.locator("#id_reference")
    reference_input.fill("not-a-doi")
    reference_input.blur()

    # Should show validation error
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("DOI")

    # Enter valid DOI
    reference_input.fill("10.1234/test.doi")
    reference_input.blur()

    # Should now show ready message
    expect(validation_msg).to_contain_text("ready to submit")


def test_submit_button_enabled_when_valid(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that submit button enables when all spectrum fields are valid."""
    page = spectrum_form_page

    submit_btn = page.locator("#submit-btn")

    # Initially disabled (no file uploaded)
    expect(submit_btn).to_be_disabled()

    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    # Still disabled after column selection (missing required fields)
    expect(submit_btn).to_be_disabled()

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')

    # Fill spectrum fields one by one
    category_select.select_option("f")
    expect(submit_btn).to_be_disabled()  # still missing subtype, owner, source

    subtype_select.select_option("bp")
    expect(submit_btn).to_be_disabled()  # still missing owner, source

    owner_input.fill("Test Filter")
    expect(submit_btn).to_be_disabled()  # still missing source

    page.locator("#id_source").fill("Test source")

    # Submit button is enabled once all spectrum fields and source are valid
    # (confirmation checkbox is validated at HTML form submission, not by JS)
    expect(submit_btn).to_be_enabled()

    # Validation message should show ready state
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("ready to submit")


# --- Spectrum Card Actions ---


def test_remove_spectrum(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test removing a spectrum card."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1, 2])

    # Should have 2 cards
    expect(page.locator(".spectrum-card")).to_have_count(2)

    # Remove first card
    first_card = page.locator(".spectrum-card").first
    remove_btn = first_card.locator('[id^="remove-btn-"]')
    remove_btn.click()

    # Should have 1 card left
    expect(page.locator(".spectrum-card")).to_have_count(1)


def test_remove_all_spectra_hides_source_fields(
    spectrum_form_page: Page, sample_csv_file: Path
) -> None:
    """Test that removing all spectra hides the source fields."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    # Source fields should be visible
    expect(page.locator("#global-source-fields")).to_be_visible()

    # Remove the only card
    remove_btn = page.locator('[id^="remove-btn-"]').first
    remove_btn.click()

    # Source fields should be hidden
    expect(page.locator("#global-source-fields")).not_to_be_visible()


# --- Chart Interaction Tests ---


def test_chart_click_sets_peak(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that clicking on the chart sets the peak marker."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    peak_badge = card.locator('[id^="peak-badge-"]')

    # Select normalizable category
    category_select.select_option("p")
    subtype_select.select_option("ex")

    # Peak badge should show auto-detected peak
    expect(peak_badge).to_contain_text("Peak:")

    # Chart should be visible
    chart = card.locator(".highcharts-container")
    expect(chart).to_be_visible()

    # Click on the chart area (center)
    chart.click()

    # Peak badge should update (value may change based on click location)
    expect(peak_badge).to_contain_text("nm")


# --- Form Submission Test ---


def test_form_ready_for_submission(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that form shows ready state when all fields are valid.

    Note: Actual form submission is tested separately. This test verifies
    that the client-side validation and JSON data generation work correctly.
    """
    page = spectrum_form_page

    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[3])  # Filter column

    card = page.locator(".spectrum-card").first
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')

    # Fill all required fields
    category_select.select_option("f")
    subtype_select.select_option("bp")
    owner_input.fill("Test BP Filter 525")
    page.locator("#id_source").fill("E2E Test")

    # Verify validation message shows ready
    validation_msg = page.locator("#validation-message")
    expect(validation_msg).to_contain_text("ready to submit")

    # Verify hidden JSON field is populated
    spectra_json = page.locator("#id_spectra_json")
    json_value = spectra_json.input_value()
    assert json_value != "[]", "spectra_json should contain data"
    assert "Test BP Filter 525" in json_value, "JSON should contain owner name"
    assert '"category":"f"' in json_value, "JSON should contain category"
    assert '"subtype":"bp"' in json_value, "JSON should contain subtype"

    # Submit button should be enabled
    submit_btn = page.locator("#submit-btn")
    expect(submit_btn).to_be_enabled()


# --- Status Indicator Tests ---


def test_status_indicators_update(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that status indicators (✓/!) update as fields are filled."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first
    status_icon = card.locator('[id^="status-icon-"]')

    # Initially should show warning
    expect(status_icon).to_contain_text("!")

    # Fill required fields
    category_select = card.locator('[id^="category-select-"]')
    subtype_select = card.locator('[id^="subtype-select-"]')
    owner_input = card.locator('[id^="owner-input-"]')

    category_select.select_option("f")
    subtype_select.select_option("bp")
    owner_input.fill("Test Filter")

    # Should now show checkmark
    expect(status_icon).to_contain_text("✓")


# --- Label Bold State Tests ---


def test_labels_bold_when_empty(spectrum_form_page: Page, sample_csv_file: Path) -> None:
    """Test that required field labels are bold when empty."""
    page = spectrum_form_page
    _upload_csv_file(page, sample_csv_file)
    _select_columns(page, wavelength_col=0, data_cols=[1])

    card = page.locator(".spectrum-card").first

    # Category label should be bold initially
    category_label = card.locator('[id^="category-label-"]')
    expect(category_label).to_have_css("font-weight", "700")

    # Select a category
    category_select = card.locator('[id^="category-select-"]')
    category_select.select_option("f")

    # Category label should no longer be bold
    expect(category_label).to_have_css("font-weight", "400")

    # Subtype label should be bold (not yet selected for multi-option categories)
    subtype_label = card.locator('[id^="subtype-label-"]')
    expect(subtype_label).to_have_css("font-weight", "700")


# --- European CSV Format Tests ---


@pytest.fixture
def european_csv_file() -> Iterator[Path]:
    """Create a temporary European-format CSV file (semicolon delimiter, comma decimal)."""
    content = """Wavelength;Relative intensity
476,0599976;11,14449883
476,9599915;5,403448582
478,0299988;2,37003231
479,1699982;1,089500189
480,1300049;0,500396013
481,2699966;0,294151127
482,2299957;0,173027024
483,3699951;0,102012008
484,3300018;0,077506006
485,4699936;0,06600528
486,4300003;0,054604368
487,5699921;0,046903752
488,5300064;0,044503561
489,6699982;0,048103846
490,6300049;0,062004961
491,6199951;0,098407871
492,7699966;0,202816248
493,7299957;0,519941568
494,8699951;1,0
496,0100021;0,873469949
497,1500015;0,461537003
498,1100006;0,293451279
499,25;0,200416033
500,2099991;0,143111453
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        path = Path(f.name)
    yield path
    path.unlink(missing_ok=True)


def test_european_csv_format_parsing(spectrum_form_page: Page, european_csv_file: Path) -> None:
    """Test that European CSV format (semicolon delimiter, comma decimal) is parsed correctly."""
    page = spectrum_form_page
    _upload_csv_file(page, european_csv_file)

    # Column picker should be visible
    column_picker = page.locator("#column-picker-container")
    expect(column_picker).to_be_visible()

    # Should show both column headers
    headers = column_picker.locator("th.column-header")
    expect(headers).to_have_count(2)  # Wavelength + Relative intensity

    # Check that numeric values are displayed in the table
    # (not "-" which indicates parsing failure)
    # The table should show numeric values with 4 decimal places
    first_cell = column_picker.locator("table tbody tr").first.locator("td").first
    cell_text = first_cell.text_content()

    # Should be a number around 476.06, not "-" (which indicates NaN)
    assert cell_text != "-", "First cell should contain a numeric value, not '-'"
    assert "476." in cell_text, f"First cell should contain wavelength ~476, got: {cell_text}"

    # Verify we can select columns and continue
    _select_columns(page, wavelength_col=0, data_cols=[1])

    # Spectrum card should appear
    expect(page.locator(".spectrum-card")).to_have_count(1)
