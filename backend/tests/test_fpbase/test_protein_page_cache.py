"""Invalidating a cached protein page must reach every visitor's copy of it."""

from __future__ import annotations

import pytest
from django.core.cache import cache
from django.test import Client, RequestFactory

from fpbase.util import uncache_protein_page
from proteins.factories import ProteinFactory

pytestmark = pytest.mark.django_db


@pytest.fixture(autouse=True)
def _clear_cache():
    cache.clear()
    yield
    cache.clear()


def test_uncache_protein_page_clears_every_cookie_variant():
    protein = ProteinFactory(name="OldName")
    url = protein.get_absolute_url()
    editor, visitor = Client(), Client()
    editor.cookies["sessionid"] = "editor"
    visitor.cookies["sessionid"] = "visitor"
    for client in (editor, visitor):
        client.get(url)  # (sets the csrftoken cookie, which changes the cache variant)
        assert b"OldName" in client.get(url).content

    type(protein).objects.filter(id=protein.id).update(name="NewName")
    assert b"OldName" in visitor.get(url).content  # (still cached)

    cookie = "; ".join(f"{k}={v.value}" for k, v in editor.cookies.items())
    request = RequestFactory().get(url, HTTP_COOKIE=cookie)
    uncache_protein_page(protein.slug, request)
    for client in (editor, visitor):
        assert b"NewName" in client.get(url).content


def test_uncache_protein_page_leaves_other_proteins_cached():
    protein, other = ProteinFactory(name="OldName"), ProteinFactory(name="OtherOld")
    client = Client()
    for _ in range(2):  # (the first response sets the csrftoken cookie)
        client.get(other.get_absolute_url())
    type(other).objects.filter(id=other.id).update(name="OtherNew")

    uncache_protein_page(protein.slug, RequestFactory().get("/"))
    assert b"OtherOld" in client.get(other.get_absolute_url()).content
