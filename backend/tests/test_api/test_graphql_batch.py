"""The GraphQL endpoint takes one operation per request."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.django_db

QUERY = {"query": "{ proteins { id } }"}


def test_single_operation(client):
    response = client.post("/graphql/", QUERY, content_type="application/json")
    assert response.status_code == 200


def test_list_of_operations_is_rejected(client):
    response = client.post("/graphql/", [QUERY, QUERY], content_type="application/json")
    assert response.status_code == 400
    assert response.json()["errors"]


def test_batch_route_is_not_mounted(client):
    response = client.post("/graphql/batch/", [QUERY, QUERY], content_type="application/json")
    assert response.status_code == 404
