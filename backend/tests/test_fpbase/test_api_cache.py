"""API responses are cached until the data changes."""

from __future__ import annotations

import json
from contextlib import contextmanager

import pytest
import reversion
from django.core.cache import cache
from django.db import connection, transaction
from django.test.utils import CaptureQueriesContext
from reversion.models import Version

from fpbase import cache_utils, views
from fpbase.cache_utils import get_data_version, get_model_version
from proteins.factories import ProteinFactory, StateFactory
from proteins.models import FluorescenceMeasurement, Protein

pytestmark = pytest.mark.django_db


@contextmanager
def _assert_num_selects(num: int):
    # (rather than counting every query: with ATOMIC_REQUESTS, as in CI, each request
    # also runs SAVEPOINT and RELEASE SAVEPOINT)
    with CaptureQueriesContext(connection) as ctx:
        yield
    selects = [q["sql"] for q in ctx.captured_queries if q["sql"].startswith("SELECT")]
    assert len(selects) == num, selects


def _graphql(client, query: str, **variables):
    body = json.dumps({"query": query, "variables": variables})
    return client.post("/graphql/", body, content_type="application/json")


def test_rest_list_is_cached_until_data_changes(client):
    protein = ProteinFactory(name="OldName")
    url = "/api/proteins/?format=json"

    response = client.get(url)
    assert [p["name"] for p in response.json()] == ["OldName"]
    # clients and the CDN can't see the data version: they get a short lifetime
    assert response["Cache-Control"] == "max-age=600"

    with _assert_num_selects(0):
        cached = client.get(url)
    assert cached.json() == response.json()
    assert cached["Cache-Control"] == "max-age=600"
    assert "Age" not in cached

    protein.name = "NewName"
    protein.save()
    assert [p["name"] for p in client.get(url).json()] == ["NewName"]

    # a write that sends no signal is only noticed when the cache entry expires
    Protein.objects.filter(id=protein.id).update(name="Unsignaled")
    assert [p["name"] for p in client.get(url).json()] == ["NewName"]


def test_related_models_change_the_data_version():
    protein = ProteinFactory()
    version = get_data_version()
    assert get_data_version() == version

    state = StateFactory(protein=protein)
    assert (created := get_data_version()) != version
    state.delete()
    assert get_data_version() != created


def test_measurement_edit_changes_the_data_version(client):
    state = StateFactory(ex_max=488)
    url = f"/api/proteins/{state.protein.slug}/"
    assert client.get(url).json()["states"][0]["ex_max"] == 488

    # loaded on its own (as in the admin), its `state` is a FluorState, not a State
    measurement = FluorescenceMeasurement.objects.get(state=state)
    measurement.ex_max = 500
    measurement.save()
    assert client.get(url).json()["states"][0]["ex_max"] == 500


@pytest.fixture
def real_on_commit(monkeypatch):
    """Undo the conftest fixture that invalidates without waiting for a commit."""
    monkeypatch.setattr(cache_utils, "_after_commit", transaction.on_commit)


@pytest.mark.usefixtures("real_on_commit")
def test_data_version_changes_on_commit_only(django_capture_on_commit_callbacks):
    protein = ProteinFactory()
    version = get_data_version()

    with (
        django_capture_on_commit_callbacks(execute=True) as callbacks,
        pytest.raises(RuntimeError),
        transaction.atomic(),
    ):
        protein.save()
        raise RuntimeError("rolled back")
    assert not callbacks
    assert get_data_version() == version

    with django_capture_on_commit_callbacks(execute=True):
        protein.save()
        # other requests still read the old row: nothing to invalidate yet
        assert get_data_version() == version
    assert get_data_version() != version


@pytest.mark.usefixtures("real_on_commit")
def test_viewing_a_protein_version_keeps_the_data_version(
    client, django_capture_on_commit_callbacks
):
    with reversion.create_revision():
        protein = ProteinFactory(blurb="Original blurb")
    with reversion.create_revision():
        protein.blurb = "Modified blurb"
        protein.save()
    first = Version.objects.get_for_object(protein).last()
    version = get_data_version()

    # the page reverts the revision to render it, and rolls that back
    with django_capture_on_commit_callbacks(execute=True):
        response = client.get(f"/protein/{protein.slug}/ver/{first.id}")
    assert response.status_code == 200
    assert b"Original blurb" in response.content
    assert get_data_version() == version


def test_model_versions_do_not_expire():
    version = get_model_version(Protein)
    # (with the default timeout, the version changed every 5 minutes)
    key = cache.make_key(cache_utils._model_cache_key(Protein))
    assert cache._expire_info[key] is None  # (LocMemCache)
    assert get_model_version(Protein) == version


def test_graphql_response_is_cached_until_data_changes(client):
    protein = ProteinFactory(name="OldName")
    query = "query getProtein($id: String!) { protein(id: $id) { name } }"

    response = _graphql(client, query, id=protein.uuid)
    assert response.json() == {"data": {"protein": {"name": "OldName"}}}
    with _assert_num_selects(0):
        assert _graphql(client, query, id=protein.uuid).json() == response.json()
        # the same query as a GET shares the entry
        get = client.get("/graphql/", {"query": query, "variables": f'{{"id": "{protein.uuid}"}}'})
        assert get.json() == response.json()

    # other variables are another entry
    other = ProteinFactory(name="Other")
    assert _graphql(client, query, id=other.uuid).json()["data"]["protein"]["name"] == "Other"

    protein.name = "NewName"
    protein.save()
    assert _graphql(client, query, id=protein.uuid).json()["data"]["protein"]["name"] == "NewName"


def test_graphql_errors_and_large_responses_are_not_cached(client, monkeypatch):
    ProteinFactory(name="Appears")
    missing = '{ dye(name: "no such dye") { name } }'
    for _ in range(2):  # (a cached response would run no query)
        with _assert_num_selects(1):
            response = _graphql(client, missing)
        assert response.status_code == 200
        assert "errors" in response.json()

    # pretty-printed: "data" sorts before "errors"
    for _ in range(2):
        with _assert_num_selects(1):
            response = client.get("/graphql/", {"query": missing, "pretty": "1"})
        assert response.content.lstrip(b"{ \n").startswith(b'"data"')
        assert "errors" in response.json()

    monkeypatch.setattr(views, "GRAPHQL_CACHE_MAX_SIZE", 10)
    for _ in range(2):
        with _assert_num_selects(1):
            response = _graphql(client, "{ proteins { name } }")
        assert response.json() == {"data": {"proteins": [{"name": "Appears"}]}}
