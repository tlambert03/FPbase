"""Responses are cached at the CDN, and purged there when the data changes."""

from __future__ import annotations

import json
from io import StringIO
from unittest import mock

import pytest
from django.core.cache import cache
from django.core.management import call_command

from fpbase import edge_cache
from fpbase.tasks import purge_edge_cache
from proteins.factories import ProteinFactory
from proteins.models import Protein

pytestmark = pytest.mark.django_db


@pytest.fixture
def enabled(settings, apply_async):
    # (`apply_async` is mocked: an eager task would call Cloudflare for real)
    settings.CLOUDFLARE_ZONE_ID = "zone123"
    settings.CLOUDFLARE_PURGE_TOKEN = "secret"
    settings.CANONICAL_URL = "https://www.example.org"


@pytest.fixture
def apply_async():
    with mock.patch.object(purge_edge_cache, "apply_async") as patched:
        yield patched


def _graphql_get(client, query, **params):
    return client.get("/graphql/", {"query": query, **params})


def test_off_until_configured(client, apply_async):
    assert not edge_cache.is_enabled()
    ProteinFactory()
    apply_async.assert_not_called()
    assert client.get("/api/proteins/?format=json")["Cache-Control"] == "max-age=600"
    assert "Cache-Control" not in _graphql_get(client, "{ proteins { name } }")


@pytest.mark.usefixtures("enabled")
def test_cdn_headers(client):
    ProteinFactory(name="Cached")
    expected = "public, max-age=60, s-maxage=86400"
    assert client.get("/api/proteins/?format=json")["Cache-Control"] == expected
    assert client.get("/api/proteins/cached/?format=json")["Cache-Control"] == expected

    response = _graphql_get(client, "{ proteins { name } }")
    assert response["Cache-Control"] == expected
    response = _graphql_get(client, "{ proteins { name } }")  # (from the server cache)
    assert response["Cache-Control"] == expected

    # a response that sets a cookie is not cached by the CDN
    assert "csrftoken" not in response.cookies
    assert "csrftoken" in _graphql_get(client, "{ nope }").cookies

    # errors are not for keeping
    response = _graphql_get(client, '{ dye(name: "no such dye") { name } }')
    assert "errors" in response.json()
    assert response["Cache-Control"] == "no-store"
    assert _graphql_get(client, "{ nope }")["Cache-Control"] == "no-store"
    # a POST cannot be cached by the CDN
    body = json.dumps({"query": "{ proteins { name } }"})
    assert "Cache-Control" not in client.post("/graphql/", body, content_type="application/json")


@pytest.mark.usefixtures("enabled")
def test_changes_are_purged_once_per_burst(apply_async):
    protein = ProteinFactory()
    apply_async.assert_called_once_with(countdown=edge_cache.PURGE_DELAY)

    for _ in range(3):  # (within PURGE_DELAY of the first)
        protein.save()
    assert apply_async.call_count == 1

    cache.delete(edge_cache.PURGE_PENDING_KEY)  # (the pending purge ran)
    protein.save()
    assert apply_async.call_count == 2


@pytest.mark.usefixtures("enabled")
def test_purge_request():
    with mock.patch("fpbase.tasks.requests.post") as post:
        post.return_value.ok = True
        purge_edge_cache()
    post.assert_called_once_with(
        "https://api.cloudflare.com/client/v4/zones/zone123/purge_cache",
        json={"prefixes": ["www.example.org/api/", "www.example.org/graphql/"]},
        headers={"Authorization": "Bearer secret"},
        timeout=10,
    )


@pytest.mark.usefixtures("enabled")
def test_purge_request_with_pages():
    with mock.patch("fpbase.tasks.requests.post") as post:
        post.return_value.ok = True
        purge_edge_cache(pages=True)
    prefixes, files = (c.kwargs["json"] for c in post.call_args_list)
    assert "www.example.org/protein/" in prefixes["prefixes"]
    assert "www.example.org/api/" in prefixes["prefixes"]
    # (not a prefix purge of the whole host, which would drop the static files too)
    assert files == {"files": ["https://www.example.org/"]}


@pytest.mark.usefixtures("enabled")
def test_invalidate_api_cache_command(client, apply_async):
    ProteinFactory(name="Stale")
    url = "/api/proteins/?format=json&fields=name"
    assert client.get(url).json() == [{"name": "Stale"}]
    Protein.objects.update(name="Fresh")  # (no signal: the cache cannot know)
    assert client.get(url).json() == [{"name": "Stale"}]
    apply_async.reset_mock()
    cache.delete(edge_cache.PURGE_PENDING_KEY)  # (the purge for the factory's save ran)

    out = StringIO()
    call_command("invalidate_api_cache", stdout=out)
    assert "purges are queued" in out.getvalue()
    # the API purge, and a later one that includes the pages
    apply_async.assert_called_with(kwargs={"pages": True}, countdown=120)
    assert client.get(url).json() == [{"name": "Fresh"}]
