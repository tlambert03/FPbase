"""GraphQL list queries: relations come from joins/prefetches, not a query per row."""

from __future__ import annotations

import pytest
from django.db import connection
from django.test.utils import CaptureQueriesContext

from proteins.factories import (
    MicroscopeFactory,
    OpticalConfigWithFiltersFactory,
    ProteinFactory,
    StateFactory,
)
from proteins.models import (
    Filter,
    FilterPlacement,
    OSERMeasurement,
    Protein,
    StateTransition,
)
from references.factories import AuthorFactory
from references.models import ReferenceAuthor

pytestmark = pytest.mark.django_db


def _query(client, query: str) -> tuple[dict, int]:
    with CaptureQueriesContext(connection) as ctx:
        response = client.post("/graphql/", {"query": query}, content_type="application/json")
    content = response.json()
    assert "errors" not in content, content["errors"]
    selects = [q for q in ctx.captured_queries if q["sql"].startswith("SELECT")]
    return content["data"], len(selects)


def _add_protein(i: int) -> None:
    """A protein with every relation the queries below ask for."""
    protein = ProteinFactory(name=f"Protein{i}")
    dark = StateFactory(protein=protein, name="dark")
    default = protein.states.exclude(id=dark.id).get()
    protein.default_state = default
    protein.save()
    StateTransition.objects.create(protein=protein, from_state=dark, to_state=default)
    OSERMeasurement.objects.create(
        protein=protein, percent=90, reference=protein.primary_reference
    )
    author = AuthorFactory(family=f"Author{i}")
    ReferenceAuthor.objects.get_or_create(
        reference=protein.primary_reference, author=author, defaults={"author_idx": i}
    )


QUERIES = {
    "default_state": "{ proteins { name defaultState { name exMax } } }",
    "reference_authors": "{ proteins { name primaryReference { doi authors { family } } } }",
    "transitions": "{ proteins { name transitions { fromState { name } toState { name } } } }",
    "organism_oser": "{ proteins { name parentOrganism { scientificName } oser { percent } } }",
    "states": "{ states { name protein { name } spectra { id subtype } } }",
    "back_references": "{ proteins { states { protein { states { protein { name } } } } } }",
    "connection": (
        "{ allProteins(first: 50) { edges { node {"
        " name defaultState { exMax } states { name spectra { id } } } } } }"
    ),
    "references": "{ references { doi authors { family publications { doi } } } }",
    # a relation selected with nothing but a nested relation under it: the optimizer
    # must still fetch the foreign key ("cannot be both deferred and traversed")
    "only_nested": (
        "{ proteins { primaryReference { authors { family } }"
        " parentOrganism { proteins { name } } } }"
    ),
    "states_nested": "{ states { protein { primaryReference { authors { family } } } } }",
    # `oser` and `oserMeasurements` are the same relation twice: two prefetches of it
    # with different querysets would clash
    "oser_twice": (
        "{ proteins { oser { percent } oserMeasurements { percent reference { doi } } } }"
    ),
    "oser_twice_connection": (
        "{ allProteins(first: 50) { edges { node {"
        " oser { percent } oserMeasurements { percent } } } } }"
    ),
    "spectra_reference": "{ states { spectra { id reference { doi } } } }",
    "organisms": "{ organisms { scientificName proteins { name } } }",
}


@pytest.mark.parametrize("name", QUERIES)
def test_query_count_does_not_grow_with_rows(client, name):
    _add_protein(0)
    _add_protein(1)
    _, few = _query(client, QUERIES[name])
    for i in range(2, 6):
        _add_protein(i)
    _, many = _query(client, QUERIES[name])
    assert many == few


def test_relations_in_list_queries(client):
    _add_protein(0)
    protein = Protein.objects.get(name="Protein0")
    dark, default = protein.states.get(name="dark"), protein.default_state

    data, _ = _query(client, "{ proteins { defaultState { id } } }")
    assert data["proteins"] == [{"defaultState": {"id": str(default.id)}}]
    data, _ = _query(client, "{ proteins { transitions { fromState { id } toState { id } } } }")
    assert data["proteins"][0]["transitions"] == [
        {"fromState": {"id": str(dark.id)}, "toState": {"id": str(default.id)}}
    ]
    data, _ = _query(client, QUERIES["reference_authors"])
    assert {"family": "Author0"} in data["proteins"][0]["primaryReference"]["authors"]
    # these two always failed, with "Invalid field name(s) given in select_related"
    data, _ = _query(client, QUERIES["references"])
    assert "Author0" in [a["family"] for ref in data["references"] for a in ref["authors"]]
    data, _ = _query(client, QUERIES["organisms"])
    assert [{"name": "Protein0"}] in [o["proteins"] for o in data["organisms"]]
    data, _ = _query(client, QUERIES["oser_twice"])
    assert data["proteins"][0]["oser"] == [{"percent": 90.0}]
    assert data["proteins"][0]["oserMeasurements"][0]["percent"] == 90.0


def test_missing_protein_is_null(client):
    data, _ = _query(client, '{ protein(id: "NOPE1") { name } }')
    assert data == {"protein": None}


def test_optical_config_filter_without_spectrum(client):
    config = OpticalConfigWithFiltersFactory(microscope=MicroscopeFactory())
    bare = Filter.objects.create(name="No spectrum")
    FilterPlacement.objects.create(filter=bare, config=config, path=FilterPlacement.EM)

    # (this failed with "Filter has no spectrum.")
    data, _ = _query(client, "{ opticalConfigs { filters { id spectrumId spectrum { id } } } }")
    filters = {f["id"]: f for f in data["opticalConfigs"][0]["filters"]}
    assert len(filters) == 4
    for filter_id, placement in filters.items():
        if filter_id == str(bare.id):
            assert placement == {"id": filter_id, "spectrumId": None, "spectrum": None}
        else:
            assert placement["spectrumId"] == placement["spectrum"]["id"]
