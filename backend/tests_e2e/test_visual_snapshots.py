"""Comprehensive visual snapshot tests for FPbase.

This module provides broad coverage of all main pages and their modals.
Tests are organized by page type for maintainability.

To update snapshots:
    uv run pytest backend/tests_e2e/ --assert-snapshots -n 4 --update-snapshots

To verify snapshots:
    uv run pytest backend/tests_e2e/ --assert-snapshots -n 4
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import (
    FilterFactory,
    MicroscopeFactory,
    OpticalConfigWithFiltersFactory,
    ProteinFactory,
    create_egfp,
    create_rich_microscope,
)
from proteins.models import Filter

if TYPE_CHECKING:
    from collections.abc import Callable

    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


# =============================================================================
# Static Pages
# =============================================================================


def test_home_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test home page renders correctly."""
    # Create a few proteins for the recent updates
    ProteinFactory.create_batch(3)

    page.goto(f"{live_server.url}/")
    expect(page).to_have_url(f"{live_server.url}/")

    # Wait for main content
    expect(page.locator("main, #content")).to_be_visible()

    assert_snapshot(page)


def test_about_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test about page renders correctly."""
    page.goto(f"{live_server.url}{reverse('about')}")
    expect(page).to_have_url(f"{live_server.url}{reverse('about')}")

    # Wait for page content
    expect(page.get_by_role("heading", name=re.compile("About", re.IGNORECASE))).to_be_visible()

    assert_snapshot(page)


@pytest.mark.skip(reason="Lineage page requires complex protein relationship setup")
def test_lineage_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test lineage page renders correctly."""
    # TODO: Create proper protein lineage data
    url = f"{live_server.url}{reverse('proteins:lineage-list')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for the lineage container to be visible
    expect(page.locator("#lineageOuter")).to_be_visible()

    assert_snapshot(page)


# =============================================================================
# Protein Pages
# =============================================================================


def test_protein_detail_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test protein detail page renders correctly."""
    protein = create_egfp()

    url = f"{live_server.url}{reverse('proteins:protein-detail', args=(protein.slug,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for main content
    expect(page.get_by_role("heading", name=protein.name, exact=True)).to_be_visible()

    # Wait for chart to render
    if not hasattr(assert_snapshot, "NOOP"):
        page.locator(".highcharts-series").first.wait_for(state="attached")
        page.wait_for_load_state("networkidle")
        assert_snapshot(page)


@pytest.mark.skip(reason="Share modal selector needs investigation")
def test_protein_detail_spectra_modal(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test protein detail page with spectra URL builder modal."""
    protein = create_egfp()

    url = f"{live_server.url}{reverse('proteins:protein-detail', args=(protein.slug,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for chart to render
    page.locator(".highcharts-series").first.wait_for(state="attached")

    # TODO: Find correct selector for share button
    # Click "Share" button to open the modal
    share_button = page.locator('button:has-text("Share")')
    expect(share_button).to_be_visible()
    share_button.click()

    # Wait for modal to appear
    modal = page.locator(".modal-dialog")
    expect(modal).to_be_visible()

    if not hasattr(assert_snapshot, "NOOP"):
        assert_snapshot(page)


@pytest.mark.skip(reason="Protein parent relationship not available in factory")
def test_protein_detail_mutations_modal(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test protein detail page with mutations display."""
    # TODO: Create proteins with proper parent/child relationships
    seq1 = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
    ProteinFactory.create(name="GFP1", seq=seq1)
    protein2 = ProteinFactory.create(name="GFP2", seq=seq1.replace("ELDG", "ETTG"))

    url = f"{live_server.url}{reverse('proteins:protein-detail', args=(protein2.slug,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    mutations_section = page.locator("text=/Mutations/")
    expect(mutations_section).to_be_visible()

    assert_snapshot(page)


def test_protein_table(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test protein table page renders correctly."""
    # Create proteins with fixed names to ensure consistent snapshots
    for i in range(5):
        ProteinFactory.create(name=f"Protein{i:02d}")

    url = f"{live_server.url}{reverse('proteins:table')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for table to load
    expect(page.locator("table")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


# =============================================================================
# Spectra Pages
# =============================================================================


def test_spectra_viewer_empty(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test spectra viewer initial empty state."""
    url = f"{live_server.url}{reverse('proteins:spectra')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for viewer
    expect(page.locator("#spectra-viewer")).to_be_visible()

    assert_snapshot(page)


def test_spectra_url_builder(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test spectra URL builder page."""
    create_egfp()
    url = f"{live_server.url}{reverse('proteins:spectra-url-builder')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for page content to load
    expect(page.locator("body")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


# =============================================================================
# Microscope Pages
# =============================================================================


def test_microscopes_list(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test microscopes list page."""
    # Create some microscopes with fixed names
    for i in range(3):
        microscope = MicroscopeFactory(name=f"SnapshotScope{i:02d}")
        OpticalConfigWithFiltersFactory.create_batch(2, microscope=microscope)

    url = f"{live_server.url}{reverse('proteins:microscopes')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for page content
    expect(page.locator("body")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


def test_microscope_detail(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test microscope detail page."""
    microscope = create_rich_microscope()

    url = f"{live_server.url}{reverse('proteins:microscope-detail', args=(microscope.id,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for chart to render
    expect(page.locator(".svg-container svg")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


def test_microscope_detail_with_fluorophore(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test microscope detail page with fluorophore selected."""
    create_egfp()
    microscope = create_rich_microscope()

    url = f"{live_server.url}{reverse('proteins:microscope-detail', args=(microscope.id,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for select2 to be ready
    expect(page.locator("#select2-fluor-select-container")).to_be_visible()

    # Select a fluorophore
    page.locator("#select2-fluor-select-container").click()
    search_field = page.locator(".select2-search__field")
    expect(search_field).to_be_visible()
    search_field.type("EGFP")
    highlighted = page.locator(".select2-results__option--highlighted")
    expect(highlighted).to_be_visible()
    page.keyboard.press("Enter")

    # Wait for chart to update
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


# =============================================================================
# Submission Forms (Authenticated)
# =============================================================================


def test_protein_submit_form(auth_page: Page, live_server: LiveServer, assert_snapshot: Callable) -> None:
    """Test protein submission form."""
    url = f"{live_server.url}{reverse('proteins:submit')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Wait for page to load
    expect(auth_page.locator("body")).to_be_visible()
    auth_page.wait_for_load_state("networkidle")

    assert_snapshot(auth_page)


def test_spectra_submit_form(auth_page: Page, live_server: LiveServer, assert_snapshot: Callable) -> None:
    """Test spectra submission form."""
    # Create a protein for the owner_state field
    ProteinFactory.create(name="TestProtein")

    url = f"{live_server.url}{reverse('proteins:submit-spectra')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Wait for form to be ready
    expect(auth_page.locator("#spectrum-submit-form[data-form-ready='true']")).to_be_attached()

    assert_snapshot(auth_page)


def test_microscope_create_form(auth_page: Page, live_server: LiveServer, assert_snapshot: Callable) -> None:
    """Test microscope creation form."""
    # Create some filters for selection
    FilterFactory.create_batch(3, subtype=Filter.BP)

    url = f"{live_server.url}{reverse('proteins:newmicroscope')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Wait for form to load
    expect(auth_page.locator('input[name="name"]')).to_be_visible()

    assert_snapshot(auth_page)


# =============================================================================
# Collections & Organisms
# =============================================================================


def test_collections_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test collections list page."""
    url = f"{live_server.url}{reverse('proteins:collections')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for page content
    expect(page.locator("body")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


def test_organisms_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test organisms list page."""
    url = f"{live_server.url}{reverse('proteins:organism-list')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for page content
    expect(page.locator("body")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


# =============================================================================
# Other Pages
# =============================================================================


def test_interactive_chart(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test interactive chart page."""
    # Create proteins with fixed names for consistent snapshots
    ProteinFactory.create(
        name="ChartProt01",
        agg="m",
        default_state__ext_coeff=20000,
        default_state__qy=0.5,
        default_state__ex_max=490,
        default_state__em_max=520,
    )
    ProteinFactory.create(
        name="ChartProt02",
        agg="m",
        default_state__ext_coeff=40000,
        default_state__qy=0.6,
        default_state__ex_max=550,
        default_state__em_max=580,
    )

    url = f"{live_server.url}{reverse('proteins:ichart')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for loading to complete
    loading_gif = page.locator(".loadinggif")
    expect(loading_gif).to_be_hidden(timeout=10000)

    # Wait for chart
    expect(page.locator("#mainchart")).to_be_visible()
    page.wait_for_load_state("networkidle")

    # Skip snapshot on NOOP mode (no visual testing)
    if hasattr(assert_snapshot, "NOOP"):
        return

    assert_snapshot(page)


def test_protein_comparison(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test protein comparison page."""
    seq = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKL"
    protein1 = ProteinFactory.create(name="CompGFP1", seq=seq)
    protein2 = ProteinFactory.create(name="CompGFP2", seq=seq.replace("ELF", "ETF"))

    proteins_param = f"{protein1.slug},{protein2.slug}"
    url = f"{live_server.url}{reverse('proteins:compare', args=(proteins_param,))}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for comparison to load
    expect(page.locator("text=/Mutations/")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


def test_bleaching_page(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """Test bleaching information page."""
    url = f"{live_server.url}{reverse('bleaching')}"
    page.goto(url)
    expect(page).to_have_url(url)

    # Wait for page content
    expect(page.locator("body")).to_be_visible()
    page.wait_for_load_state("networkidle")

    assert_snapshot(page)


# =============================================================================
# Admin/Management Pages (Authenticated)
# =============================================================================


@pytest.mark.parametrize("page_name", ["problems", "problems-gaps", "problems-inconsistencies"])
def test_problems_pages(auth_page: Page, live_server: LiveServer, page_name: str, assert_snapshot: Callable) -> None:
    """Test problems dashboard pages."""
    url = f"{live_server.url}{reverse(f'proteins:{page_name}')}"
    auth_page.goto(url)
    expect(auth_page).to_have_url(url)

    # Wait for page content
    expect(auth_page.locator("body")).to_be_visible()
    auth_page.wait_for_load_state("networkidle")

    assert_snapshot(auth_page)
