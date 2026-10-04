"""GraphQL queries may nest objects up to `GRAPHQL_MAX_DEPTH` levels."""

from __future__ import annotations

import pytest
from graphql import get_introspection_query

from fpbase.views import GRAPHQL_MAX_DEPTH

pytestmark = pytest.mark.django_db


def _nested(depth: int) -> str:
    """`{ proteins { states { protein { states { ... id } } } } }`, `depth` levels deep."""
    fields = ["proteins", *(["states", "protein"] * depth)][:depth]
    query = "id"
    for field in reversed(fields):
        query = f"{field} {{ {query} }}"
    return f"{{ {query} }}"


def _post(client, query: str):
    return client.post("/graphql/", {"query": query}, content_type="application/json")


def test_query_at_max_depth_runs(client):
    response = _post(client, _nested(GRAPHQL_MAX_DEPTH))
    assert response.status_code == 200
    assert "errors" not in response.json()


def test_query_beyond_max_depth_is_rejected(client):
    response = _post(client, _nested(GRAPHQL_MAX_DEPTH + 1))
    assert response.status_code == 400
    (error,) = response.json()["errors"]
    assert f"exceeds maximum operation depth of {GRAPHQL_MAX_DEPTH}" in error["message"]


def test_introspection_is_not_depth_limited(client):
    response = _post(client, get_introspection_query())
    assert response.status_code == 200
    assert "errors" not in response.json()


def test_standard_validation_still_applies(client):
    response = _post(client, "{ proteins { notAField } }")
    assert response.status_code == 400
    assert "Cannot query field 'notAField'" in response.json()["errors"][0]["message"]


def test_unknown_fragment_is_rejected(client):
    response = _post(client, "{ spectrum(id: 1) { owner { id ...FluorophoreParts } } }")
    assert response.status_code == 400
    assert "Unknown fragment 'FluorophoreParts'" in response.json()["errors"][0]["message"]


def test_fragment_cycle_is_rejected(client):
    query = """
        { proteins { ...A } }
        fragment A on Protein { id ...B }
        fragment B on Protein { id ...A }
    """
    response = _post(client, query)
    assert response.status_code == 400
    assert "Cannot spread fragment 'A'" in response.json()["errors"][0]["message"]
