"""A protein with status "hidden" is not shown to the public, on pages or in the APIs."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest
from django.core.cache import cache
from django.urls import reverse

from proteins.factories import ProteinFactory
from proteins.models import Lineage, OSERMeasurement, Protein, ProteinCollection
from references.factories import AuthorFactory
from references.models import ReferenceAuthor
from tests.test_users.factories import UserFactory

if TYPE_CHECKING:
    from django.test import Client

pytestmark = pytest.mark.django_db

NAME, ALIAS = "HiddenCanaryFP", "CanaryAlias"


@pytest.fixture(autouse=True)
def _clear_cache():
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def hidden() -> Protein:
    """A hidden protein with states, spectra, an organism, a reference, a lineage, ..."""
    parent = Lineage.objects.create(protein=ProteinFactory(name="VisibleParentFP"))
    protein = ProteinFactory(name=NAME, aliases=[ALIAS], status=Protein.STATUS.hidden)
    Lineage.objects.create(protein=protein, parent=parent, mutation="A2G")
    collection = ProteinCollection.objects.create(name="public collection", owner=UserFactory())
    collection.proteins.add(protein, parent.protein)
    reference = protein.primary_reference
    OSERMeasurement.objects.create(protein=protein, reference=reference, percent=50)
    protein.default_state.spectra.update(reference=reference)
    return protein


def _leaks(response) -> list[str]:
    if response.status_code == 404:  # (the "not found" page repeats the requested url)
        return []
    content = b"".join(response.streaming_content) if response.streaming else response.content
    text = (content.decode(errors="replace") + response.headers.get("Location", "")).lower()
    return [marker for marker in (NAME, ALIAS) if marker.lower() in text]


def _api_urls(protein: Protein) -> dict[str, str]:
    spectrum = protein.default_state.spectra.first()
    return {
        "proteins": "/api/proteins/?format=json",
        "proteins by name": f"/api/proteins/?name__icontains={NAME[:12]}&format=json",
        "proteins by status": "/api/proteins/?status=hidden&format=json",
        "protein": f"/api/proteins/{protein.slug}/?format=json",
        "protein by alias": f"/api/proteins/{ALIAS}/?format=json",
        "basic": "/api/proteins/basic/?format=json",
        "states": "/api/proteins/states/?format=json",
        "spectra": "/api/proteins/spectra/?format=json&limit=100",
        "table": reverse("api:protein-table-api"),
        "spectra-list": "/api/spectra-list/",
        "search-index": "/api/search-index/",
        "spectrum": f"/api/spectrum/{spectrum.id}/?format=json",
        "spectrum list": "/api/spectrum/?format=json",
    }


def _page_urls(protein: Protein) -> dict[str, str]:
    state = protein.default_state
    spectrum = state.spectra.first()
    collection = ProteinCollection.objects.get()
    parent = Protein.objects.get(name="VisibleParentFP")
    reference = protein.primary_reference
    author = AuthorFactory()
    ReferenceAuthor.objects.create(reference=reference, author=author, author_idx=0)
    return {
        "autocomplete protein": reverse("proteins:protein-autocomplete") + f"?q={NAME[:8]}",
        "autocomplete state": reverse("proteins:state-autocomplete") + f"?q={NAME[:8]}",
        "autocomplete lineage": reverse("proteins:lineage-autocomplete") + f"?q={NAME[:8]}",
        "detail": protein.get_absolute_url(),
        "detail by alias": f"/protein/{ALIAS}/",
        "detail by id": f"/protein/{protein.uuid}/",
        "history": reverse("proteins:protein-history", args=(protein.slug,)),
        "bleach form": reverse("proteins:protein-bleach-form", args=(protein.slug,)),
        "spectra image": reverse("proteins:spectra-img", args=(protein.slug, "svg")),
        "spectra csv": reverse("proteins:spectra_csv") + f"?q={spectrum.id}",
        "spectra of owner": f"/spectra/{state.slug}",
        "widget": reverse("proteins:widget-detail", args=(protein.slug,)),
        "compare": reverse("proteins:compare"),  # (the proteins of the session's list)
        "lineage json": reverse("proteins:get-lineage", args=(protein.slug,)),
        "lineage json (parent)": reverse("proteins:get-lineage", args=(parent.slug,)),
        "lineage json (all)": reverse("proteins:get-lineage"),
        "lineage json (organism)": reverse(
            "proteins:get-org-lineage", args=(parent.parent_organism.pk,)
        ),
        "parent detail": parent.get_absolute_url(),
        "organism": reverse("proteins:organism-detail", args=(protein.parent_organism.pk,)),
        "reference": reference.get_absolute_url(),
        "reference list": reverse("reference-list"),
        "author": reverse("references:author-detail", args=(author.pk,)),
        "collection": collection.get_absolute_url(),
        "collection json": collection.get_absolute_url() + "?format=json",
        "collection csv": collection.get_absolute_url() + "?format=csv",
        "activity": reverse("proteins:activity"),
        "search by name": reverse("proteins:search") + f"?name__icontains={NAME[:12]}",
        "fret": reverse("proteins:fret"),
        "problems gaps": reverse("proteins:problems-gaps"),
        "problems inconsistencies": reverse("proteins:problems-inconsistencies"),
    }


GRAPHQL = {
    "proteins": "{ proteins { name slug aliases } }",
    "allProteins": "{ allProteins { edges { node { name slug aliases } } } }",
    "allProteins by name": '{ allProteins(name_Icontains: "Canary") { edges { node { name } } } }',
    "protein by slug": '{ protein(slug: "hiddencanaryfp") { name slug } }',
    "protein by name": '{ protein(name: "HiddenCanaryFP") { name slug } }',
    "protein by alias": '{ protein(name: "CanaryAlias") { name slug } }',
    "protein by id": '{ protein(id: "%(uuid)s") { name slug } }',
    "states": "{ states { name slug protein { name } } }",
    "state": "{ state(id: %(state_id)s) { name slug protein { name } } }",
    "spectra owners": "{ spectra { id owner { name slug } } }",
    "spectrum": "{ spectrum(id: %(spectrum_id)s) { id owner { name slug } } }",
    "organisms": "{ organisms { proteins { name slug } } }",
    "organism": "{ organism(id: %(organism_id)s) { proteins { name slug } } }",
    "references": "{ references { proteins { edges { node { name } } } "
    "primaryProteins { edges { node { name } } } } }",
    "references oser": "{ references { oserMeasurements { protein { name slug } } } }",
    "references spectra": "{ references { spectra { owner { name slug } } } }",
}


def _graphql_params(protein: Protein) -> dict[str, object]:
    state = protein.default_state
    return {
        "uuid": protein.uuid,
        "state_id": state.id,
        "spectrum_id": state.spectra.first().id,
        "organism_id": protein.parent_organism.pk,
    }


def test_hidden_protein_is_not_in_rest_api(client: Client, hidden: Protein):
    client.force_login(UserFactory())  # (some endpoints are for logged-in users)
    leaks = {}
    for label, url in _api_urls(hidden).items():
        response = client.get(url)
        assert response.status_code < 500, (label, url)
        if found := _leaks(response):
            leaks[label] = found
    assert not leaks, f"shown by: {sorted(leaks)}"


def test_hidden_protein_is_not_on_site_pages(client: Client, hidden: Protein):
    client.force_login(UserFactory())  # (some pages are for logged-in users)
    session = client.session
    session["comparison"] = [hidden.slug, "visibleparentfp"]
    session.save()
    leaks = {}
    for label, url in _page_urls(hidden).items():
        # (as the scripts of the site's pages request them)
        response = client.get(url, headers={"x-requested-with": "XMLHttpRequest"})
        assert response.status_code < 500, (label, url)
        if found := _leaks(response):
            leaks[label] = found
    assert not leaks, f"shown by: {sorted(leaks)}"


def test_hidden_protein_is_not_in_graphql(client: Client, hidden: Protein):
    leaks = {}
    for label, query in GRAPHQL.items():
        body = json.dumps({"query": query % _graphql_params(hidden) if "%(" in query else query})
        response = client.post("/graphql/", body, content_type="application/json")
        assert response.status_code == 200, (label, response.content[:300])
        if found := _leaks(response):
            leaks[label] = found
    assert not leaks, f"shown by: {sorted(leaks)}"


def test_hidden_protein_is_shown_to_its_creator_and_staff(client: Client, hidden: Protein):
    creator = UserFactory()
    Protein.objects.filter(id=hidden.id).update(created_by=creator)
    for user in (creator, UserFactory(is_staff=True)):
        client.force_login(user)
        assert NAME in client.get(hidden.get_absolute_url()).content.decode()
