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
