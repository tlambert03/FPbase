"""End-to-end tests for GraphQL ETag caching behavior."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import create_egfp

if TYPE_CHECKING:
    from pathlib import Path

    from playwright.sync_api import Page, Response
    from pytest_django.live_server_helper import LiveServer


@pytest.fixture
def spectra_url(live_server: LiveServer) -> str:
    """Get spectra viewer URL with test data created."""
    create_egfp()
    return f"{live_server.url}{reverse('proteins:spectra')}"


def test_graphql_etag_cache_workflow(persistent_page: Page, spectra_url: str, tmp_path: Path) -> None:
    """Test that GraphQL queries use ETags properly: 200 -> 304 on reload."""
    page = persistent_page

    # Track GraphQL responses
    # Register handler BEFORE navigation to catch all requests
    graphql_responses: list[Response] = []
    page.on("response", lambda r: ("/graphql/" in r.url and graphql_responses.append(r)))

    # Navigate to spectra viewer and wait for idle
    page.goto(spectra_url)
    expect(page).to_have_url(spectra_url)
    page.wait_for_load_state("networkidle")
    # Should have at least 2 GraphQL query
    n_responses: int = len(graphql_responses)
    assert n_responses >= 2, "Expected at least 2 GraphQL queries on initial load"
    assert [r.status for r in graphql_responses] == [200] * n_responses
    assert all([((etag := r.headers.get("etag")) and etag.startswith('W/"')) for r in graphql_responses])
    assert all(["operationName=" in r.url for r in graphql_responses])
    assert "_FPB_SpectraList" in graphql_responses[0].url
    assert "_FPB_OpticalConfigList" in graphql_responses[1].url

    from rich import print

    print()
    print("---------- First GraphQL Response Headers ----------")
    print(graphql_responses[0].headers)
    graphql_responses.clear()
    # Reload page to test ETag caching behavior
    page.reload()
    page.wait_for_load_state("networkidle")

    # all requests should now have '"if-none-match" headers and return 304
    print()
    print("---------- GraphQL Request Headers on Reload ----------")
    print(graphql_responses[0].request.headers)
    assert all([(r.request.headers.get("if-none-match") is not None) for r in graphql_responses])
    assert (n_responses := len(graphql_responses)) >= 2, "Expected at least 2 GraphQL queries on initial load"
    assert [r.status for r in graphql_responses] == [304] * n_responses


def test_graphql_queries_include_operation_names(page: Page, spectra_url: str) -> None:
    """Test that GraphQL queries include operationName parameter for ETag support."""
    graphql_responses: list[Response] = []

    def handle_response(response: Response) -> None:
        if "/graphql/" in response.url:
            graphql_responses.append(response)

    # Register handler BEFORE navigation
    page.on("response", handle_response)

    # Navigate to spectra viewer
    page.goto(spectra_url)
    expect(page).to_have_url(spectra_url)
    page.wait_for_load_state("networkidle")

    # Should have at least 1 GraphQL query
    assert len(graphql_responses) >= 1, "Expected at least 1 GraphQL query"

    # Verify all GraphQL queries have operationName parameter with _FPB_ prefix
    for response in graphql_responses:
        assert "operationName=" in response.url, f"GraphQL query missing operationName parameter: {response.url}"
        assert "_FPB_" in response.url, f"GraphQL query missing _FPB_ operation name prefix: {response.url}"
