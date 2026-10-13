"""Adding a reference or excerpt by DOI from the protein page."""

from __future__ import annotations

from typing import TYPE_CHECKING

import httpx
import pytest
from django.urls import reverse

from proteins.extrest import entrez
from proteins.factories import ProteinFactory
from references.models import Reference
from tests.test_users.factories import UserFactory

if TYPE_CHECKING:
    from proteins.models import Protein

pytestmark = pytest.mark.django_db

AJAX = {"HTTP_X_REQUESTED_WITH": "XMLHttpRequest"}


@pytest.fixture
def protein(client) -> Protein:
    client.force_login(UserFactory(is_staff=True))
    return ProteinFactory()


def _crossref_404(doi: str):
    request = httpx.Request("GET", f"https://api.crossref.org/works/{doi}")
    response = httpx.Response(404, request=request)
    raise httpx.HTTPStatusError("Not Found", request=request, response=response)


def test_add_reference(client, protein: Protein):
    url = reverse("proteins:add_protein_reference", args=(protein.slug,))
    response = client.post(url, {"reference_doi": " 10.1038/NMETH.2413 "}, **AJAX)
    assert response.json() == {"status": "success"}
    assert protein.references.filter(doi="10.1038/nmeth.2413").exists()


def test_add_reference_rejects_invalid_doi(client, protein: Protein, monkeypatch):
    monkeypatch.setattr(entrez, "doi_lookup", lambda doi: pytest.fail("looked up"))
    url = reverse("proteins:add_protein_reference", args=(protein.slug,))
    doi = "https://www.fpbase.org/protein/dtomato/"
    response = client.post(url, {"reference_doi": doi}, **AJAX)
    assert response.status_code == 400
    assert response.json() == {"status": "failed", "msg": "Not a valid DOI string"}


@pytest.mark.parametrize(
    ("view", "data"),
    [
        ("proteins:add_protein_reference", {"reference_doi": "10.1234/nope"}),
        (
            "proteins:add_protein_excerpt",
            {"excerpt_doi": "10.1234/nope", "excerpt_content": "bright"},
        ),
    ],
)
def test_doi_unknown_to_crossref(client, protein: Protein, monkeypatch, view, data):
    monkeypatch.setattr(entrez, "doi_lookup", _crossref_404)
    response = client.post(reverse(view, args=(protein.slug,)), data, **AJAX)
    assert response.status_code == 400
    assert response.json() == {"status": "failed", "msg": "DOI not found at Crossref"}
    assert not Reference.objects.filter(doi="10.1234/nope").exists()
