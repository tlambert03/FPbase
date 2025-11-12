"""End-to-end tests for GraphQL ETag caching behavior.

These tests verify browser behavior with ETags:
- Chrome: fetch() automatically sends If-None-Match (uses HTTP cache)
- Safari: fetch() does NOT send If-None-Match (does not use HTTP cache)

CRITICAL: We CANNOT use Playwright's network monitoring (page.on('response'))
because it disables the HTTP cache. Instead, we use Django logging to capture
what the server receives.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pytest
from django.urls import reverse

from proteins.factories import create_egfp

if TYPE_CHECKING:
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


@pytest.fixture
def spectra_url(live_server: LiveServer) -> str:
    """Get spectra viewer URL with test data created."""
    create_egfp()
    return f"{live_server.url}{reverse('proteins:spectra')}"


@pytest.fixture
def graphql_request_logger(caplog):
    """Capture Django logs for GraphQL ETag requests."""
    caplog.set_level(logging.INFO, logger="fpbase.views")
    return caplog


@pytest.mark.only_browser("chromium")
def test_chromium_etag_behavior(
    persistent_page: Page,
    spectra_url: str,
    graphql_request_logger,
    browser_name: str,
) -> None:
    """Test that Chrome's fetch() sends If-None-Match headers and gets 304 responses."""
    if browser_name != "chromium":
        pytest.skip("This test requires Chromium")

    page = persistent_page

    # First visit
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    first_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
    assert len(first_visit_logs) >= 2, "Expected at least 2 GraphQL requests on first visit"

    graphql_request_logger.clear()

    # Second visit - Chrome should send If-None-Match and get 304s
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    cache_hits = [r for r in graphql_request_logger.records if "cache hit" in r.message.lower()]

    # Chrome's native fetch() should send If-None-Match headers
    assert len(cache_hits) >= 2, (
        f"Chrome should send If-None-Match and get 304 responses, but got {len(cache_hits)} cache hits. "
        "Chrome's ETag support may be broken."
    )


@pytest.mark.only_browser("webkit")
def test_webkit_etag_behavior(
    page: Page,
    spectra_url: str,
    graphql_request_logger,
    browser_name: str,
) -> None:
    """Test that Safari's fetch() does NOT send If-None-Match headers.

    This test verifies that the localStorage workaround in client.ts is still needed.
    We override isSafari() to test native fetch() behavior without the workaround.
    """
    if browser_name != "webkit":
        pytest.skip("This test is for WebKit/Safari only")

    # Override isSafari() to force native fetch() behavior
    page.add_init_script("""
        // Prevent localStorage workaround from activating
        window.safari = undefined;
    """)

    # First visit
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    first_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
    assert len(first_visit_logs) >= 2, "Expected at least 2 GraphQL requests on first visit"

    graphql_request_logger.clear()

    # Second visit - Safari should NOT send If-None-Match
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    second_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
    cache_hits = [r for r in graphql_request_logger.records if "cache hit" in r.message.lower()]

    # Safari's native fetch() does NOT send If-None-Match headers
    assert len(cache_hits) == 0, (
        f"Safari sent If-None-Match and got {len(cache_hits)} 304 responses! "
        "Safari's native ETag support now works - you can remove the localStorage workaround in client.ts"
    )
    assert len(second_visit_logs) >= 2, (
        "Safari should have made GraphQL requests on second visit (without If-None-Match)"
    )
