from __future__ import annotations

import pytest
from django.core.cache import cache
from django.test import Client
from rest_framework.throttling import AnonRateThrottle, UserRateThrottle

from fpbase.views import RateLimitedGraphQLView


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
