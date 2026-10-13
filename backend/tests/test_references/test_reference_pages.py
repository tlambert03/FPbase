"""Reference, author, and protein pages list references without a query per item."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from django.core.cache import cache
from django.db import connection
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from proteins.factories import ProteinFactory, StateFactory
from proteins.models import BleachMeasurement, OSERMeasurement
from references.models import Author, Reference, ReferenceAuthor

if TYPE_CHECKING:
    from django.test import Client

pytestmark = pytest.mark.django_db


def _reference(n: int, authors: list[Author]) -> Reference:
    ref = Reference(doi=f"10.1234/ref{n}", title=f"Paper {n}", year=2020)
    ref.save(skipdoi=True)
    for idx, author in enumerate(authors):
        ReferenceAuthor.objects.create(reference=ref, author=author, author_idx=idx)
    return ref


def _authors() -> list[Author]:
    # created in reverse, so neither id nor name order matches the paper's author order
    names = ["Zhang", "Adams", "Miller"]
    tag = Author.objects.count()
    authors = [Author.objects.create(family=f"{n}{tag}", given=n, initials=n[0]) for n in names]
    return authors[::-1]


def _num_queries(client: Client, url: str) -> int:
    with CaptureQueriesContext(connection) as ctx:
        response = client.get(url)
    assert response.status_code == 200
    return len(ctx)


def _author_page(n_refs: int) -> str:
    authors = _authors()
    for _i in range(n_refs):
        ref = _reference(Reference.objects.count(), authors)
        ProteinFactory(primary_reference=ref)
        ProteinFactory().references.add(ref)  # a secondary protein
    return reverse("references:author-detail", args=(authors[0].id,))


def _protein_page(n_refs: int) -> str:
    protein = ProteinFactory()
    for _ in range(n_refs):
        protein.references.add(_reference(Reference.objects.count(), _authors()))
    return reverse("proteins:protein-detail", args=(protein.slug,))


def _reference_page(n_measurements: int) -> str:
    ref = _reference(Reference.objects.count(), _authors())
    for _ in range(n_measurements):
        protein = ProteinFactory()
        protein.references.add(ref)  # shown as a secondary protein
        OSERMeasurement.objects.create(protein=protein, reference=ref, percent=90)
        state = StateFactory(protein=protein)
        BleachMeasurement.objects.create(state=state, reference=ref, rate=10)
    return reverse("references:reference-detail", args=(ref.id,))


@pytest.mark.parametrize("page", [_author_page, _protein_page, _reference_page])
def test_query_count_does_not_grow_with_items(client: Client, page) -> None:
    warmup, one, three = page(1), page(1), page(3)
    _num_queries(client, warmup)  # first-request caches
    assert _num_queries(client, three) == _num_queries(client, one)


@pytest.mark.parametrize("page", [_author_page, _protein_page, _reference_page])
def test_authors_listed_in_paper_order(client: Client, page) -> None:
    html = client.get(page(1)).content.decode()
    authorlist = html[html.index("authorlist") :]
    found = [authorlist.index(n) for n in ("Miller", "Adams", "Zhang")]
    assert found == sorted(found)


def test_reference_list_query_count_does_not_grow(client: Client) -> None:
    url = reverse("reference-list")
    _reference(0, _authors())
    _num_queries(client, url)  # first-request caches
    cache.clear()  # the page is cached
    one = _num_queries(client, url)
    for n in (1, 2):
        ProteinFactory().references.add(_reference(n, _authors()))
    cache.clear()
    assert _num_queries(client, url) == one


def test_protein_secondary_reference() -> None:
    ref = _reference(0, [])
    primary, secondary, hidden = ProteinFactory.create_batch(3)
    primary.primary_reference = ref
    primary.save()
    hidden.status = "hidden"
    hidden.save()
    for protein in (primary, secondary, hidden):
        protein.references.add(ref)
    assert ref.protein_secondary_reference == [secondary]


def test_reference_page_lists_secondary_proteins(client: Client) -> None:
    html = client.get(_reference_page(1)).content.decode()
    assert "Secondary Proteins" in html
