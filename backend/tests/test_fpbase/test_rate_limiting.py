from __future__ import annotations

import pytest
from django.core.cache import cache
from django.test import Client
from rest_framework.throttling import AnonRateThrottle, UserRateThrottle

from fpbase.views import RateLimitedGraphQLView, SameOriginExemptAnonThrottle


@pytest.mark.django_db
def test_graphql_rate_limiting(monkeypatch):
    """Test that anonymous users are rate limited at 30 requests/min for GraphQL."""

    # Temporarily enable throttling for this test only by monkeypatching the view class
    monkeypatch.setattr(
        RateLimitedGraphQLView,
        "throttle_classes",
        [AnonRateThrottle, UserRateThrottle],
    )

    # Configure throttle rates (these get picked up by the throttle classes)
    monkeypatch.setitem(AnonRateThrottle.THROTTLE_RATES, "anon", "30/min")
    monkeypatch.setitem(UserRateThrottle.THROTTLE_RATES, "user", "300/min")

    # Clear cache to ensure clean slate
    cache.clear()

    client = Client()
    query = '{"query": "{ __typename }"}'
    responses = [client.post("/graphql/", query, content_type="application/json") for _ in range(35)]
    throttled = [r for r in responses if r.status_code == 429]
    assert len(throttled) >= 5, f"Expected at least 5 throttled requests, got {len(throttled)}"


@pytest.mark.django_db
def test_same_origin_exempt_throttle(monkeypatch):
    """Test that same-origin requests are exempt from rate limiting."""

    # Temporarily enable throttling for this test only
    monkeypatch.setattr(
        RateLimitedGraphQLView,
        "throttle_classes",
        [SameOriginExemptAnonThrottle, UserRateThrottle],
    )

    # Configure throttle rates
    monkeypatch.setitem(SameOriginExemptAnonThrottle.THROTTLE_RATES, "anon", "5/min")  # Low limit for testing
    monkeypatch.setitem(UserRateThrottle.THROTTLE_RATES, "user", "300/min")

    # Clear cache to ensure clean slate
    cache.clear()

    client = Client()
    query = '{"query": "{ __typename }"}'

    # Test 1: Requests WITHOUT referer should be throttled
    responses_no_referer = [client.post("/graphql/", query, content_type="application/json") for _ in range(10)]
    throttled_no_referer = [r for r in responses_no_referer if r.status_code == 429]
    assert len(throttled_no_referer) >= 5, (
        f"Expected at least 5 throttled requests without referer, got {len(throttled_no_referer)}"
    )

    # Clear cache for next test
    cache.clear()

    # Test 2: Same-origin requests (with referer) should NOT be throttled
    # In test environment, request.get_host() returns 'testserver'
    responses_with_referer = [
        client.post(
            "/graphql/",
            query,
            content_type="application/json",
            HTTP_REFERER="http://testserver/spectra/",
        )
        for _ in range(10)
    ]
    throttled_with_referer = [r for r in responses_with_referer if r.status_code == 429]
    assert len(throttled_with_referer) == 0, (
        f"Expected 0 throttled requests with same-origin referer, got {len(throttled_with_referer)}"
    )

    # Clear cache for next test
    cache.clear()

    # Test 3: Production-like URLs (with www.fpbase.org) should also work
    responses_with_prod_referer = [
        client.post(
            "/graphql/",
            query,
            content_type="application/json",
            HTTP_HOST="www.fpbase.org",
            HTTP_REFERER="https://www.fpbase.org/spectra/",
        )
        for _ in range(10)
    ]
    throttled_prod_referer = [r for r in responses_with_prod_referer if r.status_code == 429]
    assert len(throttled_prod_referer) == 0, (
        f"Expected 0 throttled requests with production referer, got {len(throttled_prod_referer)}"
    )
