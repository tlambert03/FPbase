"""End-to-end tests for FPbase using Playwright.

These tests use live_server which runs Django in a separate thread.
Database configuration and fixtures are in conftest.py.
"""

from __future__ import annotations

import os
import re
import sys
from typing import TYPE_CHECKING
from unittest.mock import patch

import pytest
from django.urls import reverse
from django_recaptcha.client import RecaptchaResponse
from playwright.sync_api import expect

from favit.models import Favorite
from proteins.factories import MicroscopeFactory, OpticalConfigWithFiltersFactory, ProteinFactory
from proteins.models import Microscope, Spectrum
from proteins.models.protein import Protein
from proteins.util.blast import _get_binary

if TYPE_CHECKING:
    from collections.abc import Callable

    from django.contrib.auth.models import AbstractUser
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer

SEQ = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
# Reverse translation of DGDVNGHKFSVSGEGEGDATYGKLTLKFICT
CDNA = "gatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcacc"


def _is_not_chromium() -> bool:
    """Check if browser specified in CLI args is not chromium.

    pytest-playwright defaults to chromium, so we only skip if a different browser
    was explicitly specified via --browser CLI option.
    """
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--browser="):
            browser = arg.split("=", 1)[1]
            return browser != "chromium"
        if arg == "--browser" and i + 1 < len(sys.argv):
            return sys.argv[i + 1] != "chromium"
    return False


def _select2_enter(selector: str, text: str, page: Page) -> None:
    """Helper to select an option in a Select2 widget by typing and selecting."""
    combo = page.locator(selector)
    combo.click()
    # Wait for search field to be ready and type
    search_field = page.locator(".select2-search__field")
    expect(search_field).to_be_visible()
    search_field.type(text)
    # Wait for and select highlighted option
    highlighted = page.locator(".select2-results__option--highlighted")
    expect(highlighted).to_be_visible()
    page.keyboard.press("Enter")


def test_main_page_loads_with_assets(live_server: LiveServer, page: Page, assert_snapshot) -> None:
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

    # Visual snapshot: capture homepage after all elements loaded
    assert_snapshot(page)


def test_spectra_viewer_loads(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test the spectra viewer page loads without console errors."""
    ProteinFactory.create()

    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    spectra_viewer = page.locator("#spectra-viewer")
    expect(spectra_viewer).to_be_attached()

    # Visual snapshot: capture spectra viewer initial state
    assert_snapshot(page)


def test_spectrum_submission_preview_manual_data(auth_page: Page, live_server: LiveServer, assert_snapshot) -> None:
    """Test spectrum submission form with manual data preview."""
    protein = ProteinFactory.create()
    protein.default_state.ex_spectrum.delete()

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)
    auth_page.locator("#id_confirmation").check()

    # Select2 autocomplete for protein
    _select2_enter("#div_id_owner_state [role='combobox']", protein.name, auth_page)

    # Switch to manual data tab and enter data
    auth_page.locator("#manual-tab").click()
    data_field = auth_page.locator("#id_data")
    expect(data_field).to_be_visible()
    data_field.fill("[[500,0.1],[505,0.5],[510,0.8],[515,0.6],[520,0.3]]")

    # Visual snapshot: form filled but before preview
    assert_snapshot(auth_page)

    # Submit for preview
    auth_page.locator('input[type="submit"]').click()

    # Wait for preview section to appear (auto-waiting)
    preview_section = auth_page.locator("#spectrum-preview-section")
    expect(preview_section).to_be_visible()

    svg = auth_page.locator("#spectrum-preview-chart svg")
    expect(svg).to_be_visible()
    expect(svg.locator("[id^='FillBetweenPolyCollection']")).to_have_count(1)

    # Visual snapshot: preview chart displayed
    if not hasattr(assert_snapshot, "NOOP"):
        auth_page.wait_for_load_state("networkidle")
        assert_snapshot(auth_page)

    # submit it!
    auth_page.get_by_text("Submit Spectrum").click()
    expect(auth_page).to_have_url(f"{live_server.url}{reverse('proteins:spectrum_submitted')}")


def test_spectrum_submission_tab_switching(auth_page: Page, live_server: LiveServer, assert_snapshot) -> None:
    """Test tab switching behavior in spectrum submission form."""
    # Create a protein so owner_state field has options
    protein = ProteinFactory.create()

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Fill out basic fields
    auth_page.locator("#id_category").select_option(Spectrum.PROTEIN)
    auth_page.locator("#id_subtype").select_option(Spectrum.EX)
    auth_page.locator("#id_confirmation").check()

    _select2_enter("#div_id_owner_state [role='combobox']", protein.name, auth_page)

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

    # Visual snapshot: manual tab active with data
    assert_snapshot(auth_page)

    # Switch back to file tab
    file_tab.click()
    expect(file_tab).to_have_class(re.compile("active"))

    # Verify submit button is present
    submit_btn = auth_page.locator('input[type="submit"]')
    expect(submit_btn).to_be_visible()

    # Visual snapshot: file tab active
    assert_snapshot(auth_page)


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


def test_microscope_create(live_server: LiveServer, auth_page: Page, assert_snapshot) -> None:
    """Test microscope creation form with optical config."""
    # Create filters that can be selected in the form
    filter_configs = OpticalConfigWithFiltersFactory.create_batch(3)
    # Extract filter names from the created configs
    ex0_name = filter_configs[0].ex_filters.first().name
    ex1_name = filter_configs[1].ex_filters.first().name
    bs_name = filter_configs[0].bs_filters.first().name
    em_name = filter_configs[0].em_filters.first().name

    # Navigate to microscopes list page
    url = f"{live_server.url}{reverse('proteins:microscopes')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Click "Create a new microscope" button
    auth_page.get_by_role("link", name="Create a new microscope").click()

    # Wait for create page to load
    create_url = f"{live_server.url}{reverse('proteins:newmicroscope')}"
    expect(auth_page).to_have_url(create_url)

    auth_page.locator('input[name="name"]').fill("test microscope")
    auth_page.locator('input[name="optical_configs-0-name"]').fill("WF Green")

    # Select first excitation filter
    # Click on the excitation filters box to open dropdown
    ex_filters_container = auth_page.locator("#div_id_optical_configs-0-ex_filters")
    ex_filters_container.click()

    # Wait for search field and type first filter name
    search_field = auth_page.locator("#div_id_optical_configs-0-ex_filters")
    expect(search_field).to_be_visible()
    search_field.type(ex0_name)
    # Wait for highlighted option to appear, then press Enter
    highlighted = auth_page.locator(".select2-results__option--highlighted")
    expect(highlighted).to_be_visible()
    auth_page.keyboard.press("Enter")

    # Verify first filter is selected
    selected_choices = auth_page.locator("#div_id_optical_configs-0-ex_filters li.select2-selection__choice")
    expect(selected_choices).to_have_count(1)

    # Wait for the dropdown results to disappear (dropdown closes after selection)
    results_dropdown = auth_page.locator(".select2-results")
    expect(results_dropdown).not_to_be_visible()

    # Click on the excitation filters box again to reopen dropdown
    ex_filters_container.click()

    # Wait for the dropdown to open by checking for results options containing our second filter
    # This ensures the dropdown is fully loaded with search results
    result_with_filter = auth_page.locator(f".select2-results__option:has-text('{ex1_name}')")
    expect(result_with_filter).to_be_visible()

    # Now type to filter the results
    search_field = auth_page.locator("#div_id_optical_configs-0-ex_filters")
    search_field.type(ex1_name)

    # Wait for the specific filter to become highlighted (ensures search completed and result is highlighted)
    highlighted_with_text = auth_page.locator(f".select2-results__option--highlighted:has-text('{ex1_name}')")
    expect(highlighted_with_text).to_be_visible()

    # Press Enter to select
    auth_page.keyboard.press("Enter")

    # Verify both filters are now selected
    selected_choices = auth_page.locator("#div_id_optical_configs-0-ex_filters li.select2-selection__choice")
    expect(selected_choices).to_have_count(2)

    # Add Dichroic Filter (bs_name)
    bs_filters_container = auth_page.locator("#div_id_optical_configs-0-bs_filters")
    bs_filters_container.click()

    # Wait for the dropdown to open with the dichroic filter
    result_with_bs = auth_page.locator(f".select2-results__option:has-text('{bs_name}')")
    expect(result_with_bs).to_be_visible()

    # Type to filter the results
    bs_search_field = auth_page.locator("#div_id_optical_configs-0-bs_filters")
    bs_search_field.type(bs_name)

    # Wait for the specific filter to become highlighted
    bs_highlighted = auth_page.locator(f".select2-results__option:has-text('{bs_name}')")
    expect(bs_highlighted).to_be_visible()
    auth_page.keyboard.press("Enter")

    # Add Emission Filter (em_name)
    em_filters_container = auth_page.locator("#div_id_optical_configs-0-em_filters")
    em_filters_container.click()

    # Wait for the dropdown to open with the emission filter
    result_with_em = auth_page.locator(f".select2-results__option:has-text('{em_name}')")
    expect(result_with_em).to_be_visible()

    # Type to filter the results
    em_search_field = auth_page.locator("#div_id_optical_configs-0-em_filters")
    em_search_field.type(em_name)

    # Wait for the specific filter to become highlighted
    em_highlighted = auth_page.locator(f".select2-results__option:has-text('{em_name}')")
    expect(em_highlighted).to_be_visible()
    auth_page.keyboard.press("Enter")

    assert_snapshot(auth_page)

    auth_page.locator('input[type="submit"]').click()
    auth_page.wait_for_load_state("networkidle")
    scope = Microscope.objects.last()
    expect(auth_page).to_have_url(f"{live_server.url}{reverse('proteins:microscope-detail', args=(scope.id,))}")
    assert_snapshot(auth_page, mask_elements=["#spectrasvg"])  # mask details of filter svgs


def test_microscope_page_with_interaction(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test microscope page with fluorophore selection and config switching."""
    protein = ProteinFactory.create()
    microscope = MicroscopeFactory(name="TestScope", id="TESTSCOPE123")
    # Only need 2 configs to test switching between them
    OpticalConfigWithFiltersFactory.create_batch(2, microscope=microscope)

    page.goto(f"{live_server.url}{reverse('proteins:microscopes')}")
    page.get_by_role("link", name=microscope.name).click()
    # Wait for the select2 widget to be ready instead of networkidle
    expect(page.locator("#select2-fluor-select-container")).to_be_visible()

    # Type fluorophore name and select
    _select2_enter("#select2-fluor-select-container", protein.name[:5], page)

    # Switch optical config - select by label (TestOC1 is the second config created)
    config_select = page.locator("#config-select")
    config_select.select_option(label="TestOC1")


def test_fret_page_loads(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test FRET page loads with donor/acceptor selection."""
    ProteinFactory(
        name="donor",
        agg="m",
        default_state__ex_max=488,
        default_state__em_max=525,
        default_state__qy=0.8,
    )
    ProteinFactory(
        name="acceptor",
        agg="m",
        default_state__ex_max=525,
        default_state__em_max=550,
        default_state__ext_coeff=55000,
        default_state__qy=0.6,
    )
    url = f"{live_server.url}{reverse('proteins:fret')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify select2 widgets are present and clickable
    _select2_enter("#select2-donor-select-container", "donor", page)
    _select2_enter("#select2-acceptor-select-container", "acceptor", page)

    # Wait for calculation to complete by checking for SVG to be rendered
    expect(page.locator("#spectra svg")).to_be_visible()

    # Verify result fields exist
    expect(page.locator("#QYD")).to_be_attached()
    expect(page.locator("#QYA")).to_be_attached()
    expect(page.locator("#overlapIntgrl")).to_be_attached()

    svg = page.locator("#spectra svg")
    expect(svg).to_be_visible()
    expect(svg.locator("g.highcharts-series")).to_have_count(5)

    # Visual snapshot: FRET calculation complete with chart
    if not hasattr(assert_snapshot, "NOOP"):
        page.wait_for_load_state("networkidle")
        assert_snapshot(page)


def test_collections_page_loads(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test collections page loads without errors."""
    url = f"{live_server.url}{reverse('proteins:collections')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Visual snapshot: collections page
    assert_snapshot(page)


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


def test_protein_table_page_loads(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test protein table page loads without errors."""
    # Create minimal data - table functionality doesn't require 10 proteins
    ProteinFactory.create_batch(3)
    url = f"{live_server.url}{reverse('proteins:table')}"
    page.goto(url)
    expect(page).to_have_url(url)


def test_interactive_chart_page(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test interactive chart page with axis selection."""
    # Create minimal data - chart interaction doesn't require 6 proteins
    ProteinFactory.create(
        name="Prot1",
        agg="d",
        default_state__ext_coeff=20000,
        default_state__qy=0.5,
        default_state__ex_max=490,
        default_state__em_max=520,
    )
    ProteinFactory.create(
        name="Prot2",
        agg="m",
        default_state__ext_coeff=40000,
        default_state__qy=0.6,
        default_state__ex_max=550,
        default_state__em_max=580,
    )
    url = f"{live_server.url}{reverse('proteins:ichart')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Visual snapshot: chart with custom axes
    if hasattr(assert_snapshot, "NOOP"):
        page.wait_for_load_state("networkidle")
        assert_snapshot(page)

    # Click X-axis radio button for quantum yield (wrapped in Bootstrap label)
    page.locator("//label[input[@id='Xqy']]").click()
    # Click Y-axis radio button for extinction coefficient (wrapped in Bootstrap label)
    page.locator("//label[input[@id='Yext_coeff']]").click()


def test_embedded_microscope_viewer(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test embedded microscope viewer with chart rendering."""
    microscope = MicroscopeFactory(name="TestScope", id="TESTSCOPE123")
    OpticalConfigWithFiltersFactory.create_batch(2, microscope=microscope)

    url = f"{live_server.url}{reverse('proteins:microscope-embed', args=(microscope.id,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Verify chart SVG rendered with content
    svg = page.locator(".svg-container svg")
    expect(svg).to_be_visible()

    # Verify the chart has rendered spectra paths
    paths = page.locator(".svg-container svg path")
    assert paths.count() > 10, f"Expected more than 10 paths, but got {paths.count()}"

    # Visual snapshot: embedded microscope viewer
    # assert_snapshot(page)


def test_protein_comparison(live_server: LiveServer, page: Page, assert_snapshot) -> None:
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

    # Visual snapshot: protein comparison with mutations
    assert_snapshot(page.get_by_text("Sequence Comparison").locator("css=+ div").screenshot())


@pytest.mark.skipif(_is_not_chromium(), reason="Timing flaky ... limiting to chrome.")
def test_advanced_search(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test advanced search with multiple filters."""
    protein = ProteinFactory.create(
        name="SearchTestGFP",
        seq=SEQ,
        default_state__ex_max=488,
        default_state__em_max=525,
    )
    protein = ProteinFactory.create(
        name="SearchTestGFP2",
        seq=SEQ + "AA",
        default_state__ex_max=600,
        default_state__em_max=650,
    )

    url = f"{live_server.url}{reverse('proteins:search')}"
    page.goto(url)
    expect(page).to_have_url(url)
    # Wait for search form to be ready
    expect(page.locator("#filter-select-0")).to_be_visible()
    assert_snapshot(page)

    # First filter: Sequence cDNA contains
    page.locator("#filter-select-0").select_option("seq")
    page.locator("#query-row-0 .operator-select").select_option("cdna_contains")
    page.locator("#id_seq__cdna_contains").fill(CDNA)

    # Add second filter row
    page.locator("#add-row-btn").click()

    # Second filter: Name starts with (row 1 doesn't have "contains" option)
    page.locator("#filter-select-1").select_option("name")
    page.locator("#query-row-1 .operator-select").select_option("istartswith")
    page.locator("#id_name__istartswith").fill(protein.name[:6])

    # Submit search
    page.locator('button[type="submit"]').first.click()
    # page.wait_for_load_state("networkidle")

    lozenges = page.locator("#ldisplay")
    expect(lozenges).to_be_visible()
    assert_snapshot(page)

    # click on table display
    page.locator("label:has(#tbutton)").click()
    table = page.locator("#tdisplay")
    expect(table).to_be_visible()

    # now verify that an exact match redirects to the protein detail page
    page.locator("#id_name__istartswith").fill(protein.name)
    page.locator('button[type="submit"]').first.click()

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


def test_blast_search(live_server: LiveServer, page: Page, assert_snapshot) -> None:
    """Test BLAST search functionality."""
    protein = ProteinFactory.create(name="BlastTestGFP", seq=SEQ)

    url = f"{live_server.url}{reverse('proteins:blast')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Fill in query sequence (with a mutation to make it interesting)
    query_input = page.locator("#queryInput")
    expect(query_input).to_be_visible()
    query_input.fill(SEQ[5:20].replace("LDG", "LG"))
    assert_snapshot(page)

    if not os.environ.get("CI"):
        try:
            _get_binary("makeblastdb")
        except Exception:
            # this is a local test without BLAST installed; skip the rest
            return

    # Submit search
    page.locator('button[type="submit"]').click()

    # Wait for results table to load
    first_result = page.locator("table tbody tr:first-child td:first-child a")
    expect(first_result).to_be_visible()

    assert_snapshot(page)

    # Verify the protein is in results
    expect(first_result).to_have_text(protein.name)

    # Click to see alignment
    first_result.click()

    # Verify we're still on the BLAST page (alignment shown in same page)
    expect(page).to_have_url(re.compile(r"/blast"))

    # Visual snapshot: BLAST alignment view
    assert_snapshot(page)


def test_favorite_button_interaction(
    auth_user: AbstractUser, auth_page: Page, live_server: LiveServer, assert_snapshot
) -> None:
    """Test favorite button interaction on protein detail page."""
    protein = Protein.objects.create(name="MyProt", seq=SEQ, uuid="XSQ4F")

    # Navigate to protein detail page
    url = f"{live_server.url}{reverse('proteins:protein-detail', args=(protein.slug,))}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Verify favorite button is present
    favorite_btn = auth_page.locator("#add_remove_favorite")

    # Verify initial state: heart icon should be hollow (.far .fa-heart)
    # Note: there are two icons (one for mobile, one for desktop) - check the visible one
    heart_icon = favorite_btn.locator("i").first
    expect(heart_icon).to_contain_class("far")
    expect(heart_icon).to_contain_class("fa-heart")

    # Verify no favorite exists in database yet
    assert Favorite.objects.get_favorite(auth_user, protein.id, "proteins.Protein") is None

    # Visual snapshot: initial state (not favorited)
    assert_snapshot(auth_page)

    # Click the favorite button
    favorite_btn.click()

    # Wait for AJAX to complete and icon to change
    # The icon should change from .far to .fas (hollow to solid)
    expect(heart_icon).to_contain_class("fas")
    expect(heart_icon).to_contain_class("fa-heart")

    # Verify backend is updated: favorite should now exist in database
    assert Favorite.objects.get_favorite(auth_user, protein.id, "proteins.Protein") is not None

    # Visual snapshot: favorited state
    assert_snapshot(auth_page)

    # Click again to unfavorite
    favorite_btn.click()

    # Should return to hollow heart
    expect(heart_icon).to_contain_class("far")

    # Verify backend is updated: favorite should be removed
    assert Favorite.objects.get_favorite(auth_user, protein.id, "proteins.Protein") is None
