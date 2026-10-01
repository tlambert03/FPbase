"""Protein detail queries should not grow with transitions or bleach measurements."""

import pytest
from django.contrib.auth.models import AnonymousUser
from django.db import connection
from django.test import RequestFactory
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory
from proteins.models import BleachMeasurement, Lineage, StateTransition
from proteins.views.protein import ProteinDetailView
from references.factories import ReferenceFactory
from references.models import Reference

pytestmark = pytest.mark.django_db


def render_detail(protein):
    request = RequestFactory().get(protein.get_absolute_url())
    request.user = AnonymousUser()
    view = ProteinDetailView()
    view.setup(request, slug=protein.slug)
    # Bypass the page cache so every measurement includes a full render.
    return view.get(request).render()


def test_detail_query_count_does_not_grow_with_measurements():
    root = Lineage.objects.create(protein=ProteinFactory(seq="", pdb=[]))
    parent = Lineage.objects.create(protein=ProteinFactory(seq="", pdb=[]), parent=root)
    protein = ProteinFactory(seq="", pdb=[])
    Lineage.objects.create(protein=protein, parent=parent)
    initial = protein.states.get()

    render_detail(protein)  # Warm shared caches such as the current Site.
    counts = []
    for index in range(4):
        state = StateFactory(protein=protein, name=f"State{index}")
        StateTransition.objects.create(
            protein=protein, from_state=initial, to_state=state, trans_wave=405
        )
        BleachMeasurement.objects.create(
            state=state, reference=protein.primary_reference, rate=10 + index
        )
        with CaptureQueriesContext(connection) as queries:
            response = render_detail(protein)
        counts.append(len(queries))
        assert response.status_code == 200
        assert f"State{index}" in response.content.decode()
        assert f"{10 + index}.0" in response.content.decode()
        assert str(parent.protein) in response.content.decode()
        assert str(root.protein) in response.content.decode()

    assert counts == [14] * len(counts), counts


def test_detail_additional_references_are_ordered_and_exclude_primary():
    protein = ProteinFactory(pdb=[])
    older = ReferenceFactory(doi="10.1038/nbt0295-151")
    newer = ReferenceFactory(doi="10.1038/373663b0")
    primary = ReferenceFactory(doi="10.1126/science.8303295")
    Reference.objects.filter(pk=older.pk).update(year=1995)
    Reference.objects.filter(pk=newer.pk).update(year=2000)
    protein.primary_reference = primary
    protein.save()
    protein.references.add(older, newer, primary)

    with CaptureQueriesContext(connection) as queries:
        response = render_detail(protein)
    assert response.context_data["additional_references"] == [newer, older]
    # The presence check and loop must share the same evaluated reference list.
    reference_queries = [q for q in queries if '"proteins_protein_references"' in q["sql"]]
    assert len(reference_queries) == 1


def test_detail_without_optional_relations():
    protein = ProteinFactory(parent_organism=None, primary_reference=None, pdb=[])
    response = render_detail(protein)
    assert response.status_code == 200
    assert response.context_data["additional_references"] == []
    assert response.context_data["has_bleach_measurements"] is False
    assert b"No photostability measurements available" in response.content
