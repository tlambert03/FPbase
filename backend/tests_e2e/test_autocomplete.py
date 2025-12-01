from __future__ import annotations

from typing import TYPE_CHECKING

from playwright.sync_api import Page, expect

if TYPE_CHECKING:
    from collections.abc import Callable

    from pytest_django.live_server_helper import LiveServer


def test_home_page_search(live_server: LiveServer, page: Page, assert_snapshot: Callable) -> None:
    """NOTE!  this actually uses the production algolia index for autocomplete suggestions."""
    page.goto(live_server.url)
    page.wait_for_load_state("networkidle")

    # The new @algolia/autocomplete-js uses different classes:
    # - .aa-Autocomplete for the wrapper
    # - .aa-Panel for the dropdown
    # - .aa-Item for each suggestion
    search_input = page.get_by_placeholder("Search")
    expect(search_input).to_be_visible()
    expect(search_input).to_be_enabled()
    assert_snapshot(page)

    # Type into the search box
    search_input.fill("egf")

    # Wait for dropdown panel to appear
    dropdown = page.locator(".aa-Panel")
    expect(dropdown).to_be_visible()
    expect(dropdown.locator(".aa-Item")).to_have_count(9, timeout=10000)
    assert_snapshot(page)
