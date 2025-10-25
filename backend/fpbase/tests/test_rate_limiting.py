from __future__ import annotations

import pytest
from django.test import Client


@pytest.mark.django_db
def test_graphql_rate_limiting():
    """Test that anonymous users are rate limited at 30 requests/min for GraphQL."""
    client = Client()
    query = '{"query": "{ __typename }"}'
    responses = [client.post("/graphql/", query, content_type="application/json") for _ in range(35)]
    throttled = [r for r in responses if r.status_code == 429]
    assert len(throttled) >= 5, f"Expected at least 5 throttled requests, got {len(throttled)}"
