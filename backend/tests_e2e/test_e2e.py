"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING
from unittest.mock import patch

import pytest
from django.urls import reverse
from django_recaptcha.client import RecaptchaResponse
from playwright.sync_api import expect

from proteins.factories import MicroscopeFactory, OpticalConfigWithFiltersFactory, ProteinFactory
from proteins.models import Spectrum
from proteins.util.blast import _get_binary

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer

SEQ = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
# Reverse translation of DGDVNGHKFSVSGEGEGDATYGKLTLKFICT
CDNA = "gatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcacc"


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


def test_spectra_viewer_loads(live_server: LiveServer, page: Page) -> None:
    """Test the spectra viewer page loads without console errors."""
    ProteinFactory.create()

    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    spectra_viewer = page.locator("#spectra-viewer")
    expect(spectra_viewer).to_be_attached()


def _select2_enter(selector: str, text: str, page: Page) -> None:
    """Helper to select an option in a Select2 widget by typing and selecting."""
    combo = page.locator(selector)
    combo.click()
    page.wait_for_selector(".select2-results", state="visible")
    page.keyboard.type(text)
    page.wait_for_selector(".select2-results__option--highlighted", state="visible")
    page.keyboard.press("Enter")
    page.wait_for_selector(".select2-results", state="hidden")


def test_spectrum_submission_preview_manual_data(auth_page: Page, live_server: LiveServer) -> None:
    """Test spectrum submission form with manual data preview."""
    protein = ProteinFactory.create()
    protein.default_state.ex_spectrum.delete()

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)

    # Select2 autocomplete for protein
    _select2_enter("#div_id_owner_state [role='combobox']", protein.name, auth_page)

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
def test_spectra_img_formats(live_server: LiveServer, page: Page, ext: str) -> None:
    """Test spectrum image generation in multiple formats."""
    protein = ProteinFactory.create()
    url = f"{live_server.url}{reverse('proteins:spectra-img', args=(protein.slug, ext))}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_spectra_img_pdf_download(live_server: LiveServer, page: Page) -> None:
    """Test spectrum PDF generation triggers download."""
    protein = ProteinFactory.create()
    url = f"{live_server.url}{reverse('proteins:spectra-img', args=(protein.slug, '.pdf'))}"

    page.goto(live_server.url)
    with page.expect_download() as download_info:
        page.evaluate(f"window.location.href = '{url}'")
    download = download_info.value
    assert download.suggested_filename.endswith(".pdf")


def test_spectra_img_with_kwargs(live_server: LiveServer, page: Page) -> None:
    """Test spectrum image generation with matplotlib kwargs."""
    protein = ProteinFactory.create()
    base_url = reverse("proteins:spectra-img", args=(protein.slug, ".svg"))
    url = f"{live_server.url}{base_url}?xlim=350,700&alpha=0.2&grid=true"
    page.goto(url)
    expect(page).to_have_url(url)


def test_microscope_page_with_interaction(live_server: LiveServer, page: Page) -> None:
    """Test microscope page with fluorophore selection and config switching."""
    protein = ProteinFactory.create()
    microscope = MicroscopeFactory(name="TestScope", id="TESTSCOPE123")
    OpticalConfigWithFiltersFactory.create_batch(4, microscope=microscope)

    page.goto(f"{live_server.url}{reverse('proteins:microscopes')}")
    page.get_by_role("link", name=microscope.name).click()
    page.wait_for_load_state("networkidle")

    # Type fluorophore name and select
    _select2_enter("#select2-fluor-select-container", protein.name[:5], page)

    # Wait for selection to complete
    page.wait_for_timeout(300)

    # Switch optical config
    config_select = page.locator("#config-select")
    config_select.click()
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    page.keyboard.press("Enter")


def test_fret_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test FRET page loads with donor/acceptor selection."""
    ProteinFactory(name="donor", agg="m", default_state__ex_max=488, default_state__em_max=525)
    ProteinFactory(
        name="acceptor",
        agg="m",
        default_state__ex_max=525,
        default_state__em_max=550,
    )
    url = f"{live_server.url}{reverse('proteins:fret')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify select2 widgets are present and clickable
    _select2_enter("#select2-donor-select-container", "donor", page)
    _select2_enter("#select2-acceptor-select-container", "acceptor", page)

    # Verify result fields exist
    expect(page.locator("#QYD")).to_be_attached()
    expect(page.locator("#QYA")).to_be_attached()
    expect(page.locator("#overlapIntgrl")).to_be_attached()

    svg = page.locator("#spectra svg")
    expect(svg).to_be_visible()
    expect(svg.locator("g.highcharts-series")).to_have_count(5)


def test_collections_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test collections page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:collections')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_problems_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test problems page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:problems')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_problems_inconsistencies_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test problems inconsistencies page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:problems-inconsistencies')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_problems_gaps_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test problems gaps page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:problems-gaps')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_protein_table_page_loads(live_server: LiveServer, page: Page) -> None:
    """Test protein table page loads without errors."""
    ProteinFactory.create_batch(10)
    url = f"{live_server.url}{reverse('proteins:table')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_interactive_chart_page(live_server: LiveServer, page: Page) -> None:
    """Test interactive chart page with axis selection."""
    ProteinFactory.create_batch(6)
    url = f"{live_server.url}{reverse('proteins:ichart')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Click X-axis radio button for quantum yield using XPath
    page.locator("//label[input[@id='Xqy']]").click()

    # Click Y-axis radio button for extinction coefficient using XPath
    page.locator("//label[input[@id='Yext_coeff']]").click()


def test_embedded_microscope_viewer(live_server: LiveServer, page: Page) -> None:
    """Test embedded microscope viewer with chart rendering."""
    microscope = MicroscopeFactory(name="TestScope", id="TESTSCOPE123")
    OpticalConfigWithFiltersFactory.create_batch(2, microscope=microscope)

    url = f"{live_server.url}{reverse('proteins:microscope-embed', args=(microscope.id,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify chart SVG rendered with content
    svg = page.locator(".svg-container svg")
    expect(svg).to_be_visible()

    paths = page.locator(".svg-container svg path")
    expect(paths).not_to_have_count(0)


def test_spectrum_submission_tab_switching(auth_page: Page, live_server: LiveServer) -> None:
    """Test tab switching behavior in spectrum submission form."""
    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Fill out basic fields
    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)

    # Wait for owner_state field to be populated and select first non-empty option
    owner_state = auth_page.locator("#id_owner_state")
    expect(owner_state).to_be_visible()
    # Wait a bit for the field to be populated after category/subtype selection
    auth_page.wait_for_timeout(500)
    # Try to select by index if options are available
    try:
        owner_state.select_option(index=1)
    except Exception:
        # If no options available, skip this step
        pass

    # Check confirmation
    auth_page.locator("#id_confirmation").check()

    # Test tab switching
    file_tab = auth_page.locator("#file-tab")
    manual_tab = auth_page.locator("#manual-tab")

    # Verify file tab is active by default
    expect(file_tab).to_have_class(re.compile("active"))

    # Switch to manual tab
    manual_tab.click()
    expect(manual_tab).to_have_class(re.compile("active"))

    # Enter manual data
    data_field = auth_page.locator("#id_data")
    expect(data_field).to_be_visible()
    expect(data_field).to_be_enabled()
    data_field.fill("[[400,0.1],[401,0.2],[402,0.3],[403,0.5],[404,0.8],[405,1.0]]")

    # Switch back to file tab
    file_tab.click()
    expect(file_tab).to_have_class(re.compile("active"))

    # Verify submit button is present
    submit_btn = auth_page.locator('input[type="submit"]')
    expect(submit_btn).to_be_visible()


def test_protein_comparison(live_server: LiveServer, page: Page) -> None:
    """Test protein comparison page shows mutations between two proteins."""
    protein1 = ProteinFactory.create(name="GFP1", seq=SEQ)
    protein2 = ProteinFactory.create(name="GFP2", seq=SEQ.replace("ELDG", "ETTG"))

    # Navigate to compare page with both proteins
    proteins_param = f"{protein1.slug},{protein2.slug}"
    url = f"{live_server.url}{reverse('proteins:compare', args=(proteins_param,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify mutations are displayed
    mutations_text = page.locator("text=/Mutations:.*L19T\\/D20T/")
    expect(mutations_text).to_be_visible()


def test_advanced_search(live_server: LiveServer, page: Page) -> None:
    """Test advanced search with multiple filters."""
    protein = ProteinFactory.create(name="SearchTestGFP", seq=SEQ)

    url = f"{live_server.url}{reverse('proteins:search')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # First filter: Sequence cDNA contains
    filter1 = page.locator("#filter-select-0")
    expect(filter1).to_be_visible()
    filter1.type("seq")
    page.keyboard.press("Tab")

    # Type to select "cDNA could contain" action
    page.keyboard.type("cdna")
    page.keyboard.press("Tab")

    # Enter cDNA value
    page.keyboard.type(CDNA)

    # Add second filter
    page.locator("#add-row-btn").click()

    # Second filter: Name contains
    filter2 = page.locator("#filter-select-1")
    expect(filter2).to_be_visible()
    filter2.type("name")
    page.keyboard.press("Tab")

    # Select "contains" action
    page.keyboard.type("cont")
    page.keyboard.press("Tab")

    # Enter partial name
    page.keyboard.type(protein.name[2:6])

    # Submit search
    page.locator('button[type="submit"]').click()

    # Should redirect to protein detail page
    expected_url = f"{live_server.url}{protein.get_absolute_url()}"
    expect(page).to_have_url(expected_url)

    # Verify protein name is shown
    heading = page.locator("h1")
    expect(heading).to_have_text(protein.name)


def test_contact_form_submission(live_server: LiveServer, page: Page) -> None:
    """Test contact form submission with reCAPTCHA.

    Verifies that the form submit button doesn't shadow the form.submit()
    method, which would break reCAPTCHA v3's JavaScript submission flow.
    """
    with patch("django_recaptcha.fields.client.submit") as mock_submit:
        mock_submit.return_value = RecaptchaResponse(is_valid=True, extra_data={"score": 0.9})

        url = f"{live_server.url}{reverse('contact')}"
        page.goto(url)
        expect(page).to_have_url(url)

        # Fill out form
        page.locator("#id_name").fill("Test User")
        page.locator("#id_email").fill("test@example.com")
        page.locator("#id_message").fill("This is a test message")

        # Wait for reCAPTCHA to load
        page.wait_for_function("typeof grecaptcha !== 'undefined'")

        # Submit form
        page.locator('input[type="submit"]').click()

        # Should redirect to thank you page
        expect(page).to_have_url(re.compile(r"/thanks/?$"))


try:
    _get_binary("makeblastdb")
    HAVE_BLAST = True
except Exception:
    HAVE_BLAST = False


@pytest.mark.skipif(not os.environ.get("CI") or not HAVE_BLAST, reason="BLAST binaries may not be installed locally")
def test_blast_search(live_server: LiveServer, page: Page) -> None:
    """Test BLAST search functionality."""
    protein = ProteinFactory.create(name="BlastTestGFP", seq=SEQ)

    url = f"{live_server.url}{reverse('proteins:blast')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Fill in query sequence (with a mutation to make it interesting)
    query_input = page.locator("#queryInput")
    expect(query_input).to_be_visible()
    query_input.fill(SEQ[5:20].replace("LDG", "LG"))

    # Submit search
    page.locator('button[type="submit"]').click()

    # Wait for results table to load
    first_result = page.locator("table tbody tr:first-child td:first-child a")
    expect(first_result).to_be_visible()

    # Verify the protein is in results
    expect(first_result).to_have_text(protein.name)

    # Click to see alignment
    first_result.click()

    # Verify we're still on the BLAST page (alignment shown in same page)
    expect(page).to_have_url(re.compile(r"/blast"))
