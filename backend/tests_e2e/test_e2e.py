"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import MicroscopeFactory, OpticalConfigWithFiltersFactory, ProteinFactory
from proteins.models import Spectrum

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer

SEQ = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"


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
    protein.default_state.ex_spectrum.delete()

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)

    # Select2 autocomplete for protein
    combo = auth_page.locator("#div_id_owner_state [role='combobox']")
    combo.click()
    auth_page.wait_for_selector(".select2-results", state="visible")
    auth_page.keyboard.type(protein.name)
    auth_page.wait_for_selector(".select2-results__option--highlighted", state="visible")
    auth_page.keyboard.press("Enter")
    auth_page.wait_for_selector(".select2-results", state="hidden")

    # Switch to manual data tab and enter data
    auth_page.locator("#manual-tab").click()
    data_field = auth_page.locator("#id_data")
    expect(data_field).to_be_visible()
    data_field.fill("[[500,0.1],[505,0.5],[510,0.8],[515,0.6],[520,0.3]]")
    auth_page.locator("#id_confirmation").check()

    # Submit for preview and wait for form to process
    submit_btn = auth_page.locator('input[type="submit"]')
    submit_btn.click()
    auth_page.wait_for_load_state("domcontentloaded")
    preview_section = auth_page.locator("#spectrum-preview-section")
    expect(preview_section).to_be_visible()

    svg = auth_page.locator("#spectrum-preview-chart svg")
    expect(svg).to_be_visible()
    expect(svg.locator("[id^='FillBetweenPolyCollection']")).to_have_count(1)

    # submit it!
    auth_page.get_by_text("Submit Spectrum").click()
    expect(auth_page).to_have_url(f"{live_server.url}{reverse('proteins:spectrum_submitted')}")


@pytest.mark.parametrize("ext", [".svg", ".png", ".jpg", ".jpeg"])
@pytest.mark.usefixtures("assert_no_console_errors")
def test_spectra_img_formats(live_server: LiveServer, page: Page, ext: str) -> None:
    """Test spectrum image generation in multiple formats."""
    protein = ProteinFactory(name="TestGFP", seq=SEQ)
    url = f"{live_server.url}{reverse('proteins:spectra-img', args=(protein.slug, ext))}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_spectra_img_pdf_download(live_server: LiveServer, page: Page) -> None:
    """Test spectrum PDF generation triggers download."""
    protein = ProteinFactory(name="TestGFP", seq=SEQ)
    url = f"{live_server.url}{reverse('proteins:spectra-img', args=(protein.slug, '.pdf'))}"

    page.goto(live_server.url)
    with page.expect_download() as download_info:
        page.evaluate(f"window.location.href = '{url}'")
    download = download_info.value
    assert download.suggested_filename.endswith(".pdf")


@pytest.mark.usefixtures("assert_no_console_errors")
def test_spectra_img_with_kwargs(live_server: LiveServer, page: Page) -> None:
    """Test spectrum image generation with matplotlib kwargs."""
    protein = ProteinFactory(name="TestGFP", seq=SEQ)
    base_url = reverse("proteins:spectra-img", args=(protein.slug, ".svg"))
    url = f"{live_server.url}{base_url}?xlim=350,700&alpha=0.2&grid=true"
    page.goto(url)
    expect(page).to_have_url(url)


@pytest.mark.usefixtures("assert_no_console_errors")
def test_microscope_page_with_interaction(live_server: LiveServer, page: Page) -> None:
    """Test microscope page with fluorophore selection and config switching."""
    protein = ProteinFactory(name="TestGFP", seq=SEQ)
    microscope = MicroscopeFactory(name="TestScope", id="TESTSCOPE123")
    OpticalConfigWithFiltersFactory.create_batch(4, microscope=microscope)

    page.goto(f"{live_server.url}{reverse('proteins:microscopes')}")
    page.get_by_role("link", name=microscope.name).click()
    page.wait_for_load_state("networkidle")

    # Wait for select2 to initialize
    fluor_combo = page.locator("#select2-fluor-select-container")
    expect(fluor_combo).to_be_visible()
    expect(fluor_combo).to_be_enabled()
    fluor_combo.click()

    # Wait for dropdown to open
    page.wait_for_selector(".select2-results", state="visible")

    # Type fluorophore name and select
    page.keyboard.type(protein.name[:5])
    page.keyboard.press("Enter")

    # Wait for selection to complete
    page.wait_for_timeout(300)

    # Switch optical config
    config_select = page.locator("#config-select")
    config_select.click()
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    page.keyboard.press("Enter")


@pytest.mark.usefixtures("assert_no_console_errors")
def test_fret_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test FRET page loads with donor/acceptor selection."""
    ProteinFactory(name="DonorFP", agg="m", default_state__ex_max=488, default_state__em_max=525)
    ProteinFactory(name="AcceptorFP", agg="m", default_state__ex_max=525, default_state__em_max=550)

    url = f"{live_server.url}{reverse('proteins:fret')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify select2 widgets are present and clickable
    donor_combo = page.locator("#select2-donor-select-container")
    acceptor_combo = page.locator("#select2-acceptor-select-container")
    expect(donor_combo).to_be_visible()
    expect(acceptor_combo).to_be_visible()

    # Verify result fields exist
    expect(page.locator("#QYD")).to_be_attached()
    expect(page.locator("#QYA")).to_be_attached()
    expect(page.locator("#overlapIntgrl")).to_be_attached()


@pytest.mark.usefixtures("assert_no_console_errors")
def test_collections_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test collections page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:collections')}"
    page.goto(url)
    expect(page).to_have_url(url)
