"""End-to-end tests for GraphQL ETag caching behavior.

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


def test_etag_behavior(
    persistent_page: Page,
    spectra_url: str,
    graphql_request_logger,
    browser_name: str,
) -> None:
    """Test that browser fetch() sends If-None-Match headers and gets 304 responses."""
    page = persistent_page

    # First visit
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    first_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
    assert len(first_visit_logs) >= 2, "Expected at least 2 GraphQL requests on first visit"

    # DEBUG: Check if localStorage has ETags after first visit
    storage_after_first = page.evaluate("""() => {
        const keys = Object.keys(localStorage);
        const etags = {};
        for (const key of keys) {
            if (key.includes('etag:')) {
                etags[key] = localStorage.getItem(key);
            }
        }
        return {
            totalKeys: keys.length,
            etagKeys: Object.keys(etags).length,
            etags: etags,
            allKeys: keys
        };
    }""")
    print(f"\n[DEBUG {browser_name}] localStorage after first visit: {storage_after_first}")

    graphql_request_logger.clear()

    # Second visit - Browser should send If-None-Match and get 304s
    page.goto(spectra_url)
    page.wait_for_load_state("networkidle")

    # DEBUG: Check if localStorage has ETags after second visit
    storage_after_second = page.evaluate("""() => {
        const keys = Object.keys(localStorage);
        const etags = {};
        for (const key of keys) {
            if (key.includes('etag:')) {
                etags[key] = localStorage.getItem(key);
            }
        }
        return {
            totalKeys: keys.length,
            etagKeys: Object.keys(etags).length,
            etags: etags,
            allKeys: keys
        };
    }""")
    print(f"\n[DEBUG {browser_name}] localStorage after second visit: {storage_after_second}")

    cache_hits = [r for r in graphql_request_logger.records if "cache hit" in r.message.lower()]

    # DEBUG: Print all second visit logs
    print(f"\n[DEBUG {browser_name}] Second visit logs:")
    for record in graphql_request_logger.records:
        print(f"  {record.levelname}: {record.message}")

    # Browser's native fetch() should send If-None-Match headers
    assert len(cache_hits) >= 2, (
        f"{browser_name} should send If-None-Match and get 304 responses, but got {len(cache_hits)} cache hits. "
        f"{browser_name}'s ETag support may be broken."
    )


# def test_webkit_etag_behavior(
#     page: Page,
#     spectra_url: str,
#     graphql_request_logger,
#     browser_name: str,
# ) -> None:
#     """Test that Safari's fetch() does NOT send If-None-Match headers.

#     This test verifies that the localStorage workaround in client.ts is still needed.
#     We override isSafari() to test native fetch() behavior without the workaround.
#     """
#     if browser_name != "webkit":
#         pytest.skip("This test is for WebKit/Safari only")

#     # Override isSafari() to force native fetch() behavior
#     page.add_init_script("""
#         // Prevent localStorage workaround from activating
#         window.safari = undefined;
#     """)

#     # First visit
#     page.goto(spectra_url)
#     page.wait_for_load_state("networkidle")

#     first_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
#     assert len(first_visit_logs) >= 2, "Expected at least 2 GraphQL requests on first visit"

#     graphql_request_logger.clear()

#     # Second visit - Safari should NOT send If-None-Match
#     page.goto(spectra_url)
#     page.wait_for_load_state("networkidle")

#     second_visit_logs = [r for r in graphql_request_logger.records if "GraphQL ETag" in r.message]
#     cache_hits = [r for r in graphql_request_logger.records if "cache hit" in r.message.lower()]

#     # Safari's native fetch() does NOT send If-None-Match headers
#     assert len(cache_hits) == 0, (
#         f"Safari sent If-None-Match and got {len(cache_hits)} 304 responses! "
#         "Safari's native ETag support now works - you can remove the localStorage workaround in client.ts"
#     )
#     assert len(second_visit_logs) >= 2, (
#         "Safari should have made GraphQL requests on second visit (without If-None-Match)"
#     )
