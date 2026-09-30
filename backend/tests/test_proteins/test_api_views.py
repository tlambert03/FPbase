"""
Tests for protein API views.

These tests ensure the API endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

import pytest
from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory
from proteins.models import Protein


class ProteinListAPIViewTests(TestCase):
    """Test the ProteinListAPIView endpoint."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests in this class."""
        for _ in range(15):
            protein = ProteinFactory()
            StateFactory(protein=protein, name="state1")

    @override_settings(DEBUG=True)
    def test_protein_list_api_query_count(self):
        """
        Test that the protein list API doesn't generate excessive queries.

        This prevents N+1 query regressions. With proper prefetch_related usage,
        the query count should be constant (4-5 queries) regardless of protein count.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.get("/api/proteins/", headers={"accept": "application/json"})

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data), 15)

        query_count = len(context.captured_queries)
        self.assertLessEqual(
            query_count,
            6,
            f"API endpoint generated {query_count} queries. "
            f"Expected ≤6 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    def test_protein_list_api_returns_expected_fields(self):
        """Test that the API returns the expected fields for each protein."""
        response = self.client.get("/api/proteins/", headers={"accept": "application/json"})

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data), 15)

        if data:
            protein = data[0]
            expected_fields = {"uuid", "name", "slug", "states", "transitions"}
            self.assertTrue(
                expected_fields.issubset(protein.keys()),
                f"Expected fields {expected_fields} not all present in {protein.keys()}",
            )

            self.assertGreater(len(protein.get("states", [])), 0, "States should be included")


class SpectraListAPIViewTests(TestCase):
    """Test the spectra-list endpoint with ETag caching."""

    @classmethod
    def setUpTestData(cls):
        """Create test proteins with states (which have spectra)."""
        cls.proteins = []
        for i in range(3):
            protein = ProteinFactory(name=f"TestProtein{i}")
            StateFactory(protein=protein, name="default")
            cls.proteins.append(protein)

    def test_spectra_list_etag_caching(self):
        """Test that ETag caching works correctly for spectra list.

        Verifies:
        1. Initial request returns 200 with ETag
        2. Second request with matching ETag returns 304
        3. After data changes, request with old ETag returns 200 with new data
        """
        # Initial request - should return 200 with ETag
        response1 = self.client.get("/api/spectra-list/")

        self.assertEqual(response1.status_code, 200)
        self.assertIn("ETag", response1)
        self.assertTrue(response1["ETag"].startswith('W/"'), "ETag should be a weak ETag")

        etag1 = response1["ETag"]
        data1 = response1.json()
        self.assertIn("data", data1)
        self.assertIn("spectra", data1["data"])
        names = {s["owner"]["name"] for s in data1["data"]["spectra"]}
        self.assertNotIn(
            "ModifiedProtein",
            names,
            "Should NOT find spectra with modified protein name",
        )

        # Second request with If-None-Match - should return 304
        response2 = self.client.get("/api/spectra-list/", headers={"If-None-Match": etag1})

        self.assertEqual(response2.status_code, 304)
        self.assertEqual(response2.content, b"")
        self.assertEqual(response2["ETag"], etag1, "304 response should include same ETag")

        # Modify a protein to invalidate cache
        protein = Protein.objects.first()
        protein.name = "ModifiedProtein"
        protein.save()

        # Third request with old ETag - should return 200 with new data
        response3 = self.client.get("/api/spectra-list/", headers={"If-None-Match": etag1})

        self.assertEqual(response3.status_code, 200, "Should return 200 after data changed")
        self.assertIn("ETag", response3)
        self.assertNotEqual(response3["ETag"], etag1, "ETag should change after data modification")

        data3 = response3.json()
        # Verify the modified protein name appears in the response
        spectra = data3["data"]["spectra"]
        # Find any spectrum owned by the modified protein
        names = {s["owner"]["name"] for s in spectra}
        self.assertIn(
            "ModifiedProtein",
            names,
            "Should find spectra with modified protein name",
        )


@pytest.mark.django_db
def test_protein_list_api_limit_offset(client):
    """limit/offset slice the (still bare) list; paging past the end gives []."""
    for _ in range(5):
        ProteinFactory()

    def get(query: str) -> list:
        return client.get(f"/api/proteins/?format=json{query}").json()

    everything = get("")
    assert len(everything) == 5
    assert get("&limit=2&offset=1") == everything[1:3]
    assert get("&limit=2&offset=1000") == []


@pytest.mark.django_db
@pytest.mark.parametrize("url", ["/api/proteins/", "/api/proteins/table-data/"])
def test_unknown_query_params_rejected(client, url):
    """Guessed params are a 400 naming the valid ones, not an unfiltered dump."""
    ProteinFactory()

    response = client.get(f"{url}?format=json&find=mCherry&per_page=2")
    assert response.status_code == 400
    error = response.json()
    assert "find, per_page" in error["detail"]
    assert "name__icontains" in error["valid_parameters"]
    assert error["docs"].endswith("/api/")

    assert client.get(f"{url}?format=json").status_code == 200


@pytest.mark.django_db
def test_protein_list_api_non_filter_params_allowed(client):
    """Search page URLs (which carry `display`) can be pasted into the API, per the docs."""
    ProteinFactory(name="KnownProtein")
    url = "/api/proteins/?format=json&name__icontains=known&display=t&limit=5&offset=0"
    response = client.get(url)
    assert response.status_code == 200
    assert [p["name"] for p in response.json()] == ["KnownProtein"]


@pytest.mark.django_db
def test_protein_list_api_name_alias(client):
    """A bare `name=` is a case-insensitive exact match."""
    ProteinFactory(name="AliasProtein")
    ProteinFactory(name="AliasProtein2")
    response = client.get("/api/proteins/?format=json&name=aliasprotein")
    assert [p["name"] for p in response.json()] == ["AliasProtein"]


@pytest.mark.django_db
def test_protein_list_api_pdb_alias(client):
    """A bare `pdb=` matches proteins with that PDB ID, in any case."""
    ProteinFactory(name="PdbProtein", pdb=["5WJ2", "2IB5"])
    ProteinFactory(name="OtherProtein", pdb=["3ADF"])
    for query in ("pdb=5WJ2", "pdb=5wj2", "pdb=2ib5,5WJ2"):
        response = client.get(f"/api/proteins/?format=json&{query}")
        assert [p["name"] for p in response.json()] == ["PdbProtein"], query


@pytest.mark.django_db
def test_spectrum_detail_api_names_protein(client):
    """A protein state's spectrum reports its protein (FPBASE-6HV)."""
    protein = ProteinFactory(name="SpecProtein")
    state = StateFactory(protein=protein, name="default")
    spectrum = state.spectra.first()
    assert spectrum is not None

    response = client.get(f"/api/spectrum/{spectrum.id}/?format=json")
    assert response.status_code == 200
    data = response.json()
    assert data["protein_name"] == "SpecProtein"
    assert data["protein_slug"] == protein.slug


@pytest.mark.django_db
def test_protein_detail_api(client, django_assert_max_num_queries):
    """A single protein by slug (in any case), in a fixed number of queries."""
    protein = ProteinFactory(name="DetailProtein")
    StateFactory(protein=protein, name="state1")
    StateFactory(protein=protein, name="state2")

    with django_assert_max_num_queries(6):
        response = client.get("/api/proteins/DetailProtein/?format=json")
    assert response.status_code == 200
    data = response.json()
    assert data["slug"] == "detailprotein"
    assert {"state1", "state2"} <= {state["name"] for state in data["states"]}

    # sibling routes are not swallowed by the slug route
    assert isinstance(client.get("/api/proteins/basic/?format=json").json(), list)
    assert client.get("/api/proteins/no-such-protein/?format=json").status_code == 404


@pytest.mark.django_db
def test_unknown_api_path_is_json_404(client):
    response = client.get("/api/no/such/endpoint/")
    assert response.status_code == 404
    assert response.json()["proteins"].endswith("/api/proteins/")

    # a missing trailing slash is a 404 too, not a redirect to one
    assert client.get("/api/no/such/endpoint").status_code == 404

    # routes declared after the api include are not shadowed
    assert client.get("/api/schema/").status_code != 404


@pytest.mark.django_db
def test_api_paths_without_trailing_slash_are_served(client):
    ProteinFactory(name="SlashProtein")

    response = client.get("/api/proteins?name=SlashProtein&format=json")
    assert response.status_code == 200
    assert [p["slug"] for p in response.json()] == ["slashprotein"]
    assert client.get("/api/proteins/slashprotein?format=json").json()["name"] == "SlashProtein"

    # a redirected POST would be resent as a GET, without the query
    response = client.post(
        "/graphql", {"query": "{ proteins { name } }"}, content_type="application/json"
    )
    assert response.status_code == 200
    assert response.json()["data"]["proteins"] == [{"name": "SlashProtein"}]

    # other pages keep the APPEND_SLASH redirect
    response = client.get("/about")
    assert response.status_code == 301
    assert response["Location"] == "/about/"


@pytest.mark.django_db
@pytest.mark.parametrize("url", ["/api/proteins/", "/api/proteins/basic/"])
def test_protein_lists_default_to_json(client, url):
    ProteinFactory(name="FormatProtein")

    # what curl and python-requests send
    response = client.get(url, headers={"accept": "*/*"})
    assert response["Content-Type"] == "application/json"
    assert response.json()[0]["name"] == "FormatProtein"

    for kwargs in ({"QUERY_STRING": "format=csv"}, {"headers": {"accept": "text/csv"}}):
        response = client.get(url, **kwargs)
        assert response["Content-Type"].startswith("text/csv")
        assert b"FormatProtein" in response.content


@pytest.mark.django_db
def test_protein_spectra_api_filters(client):
    """Filters narrow /api/proteins/spectra/, instead of being ignored (a full dump)."""
    for name in ("SpectraOne", "SpectraTwo", "SpectraThree"):
        StateFactory(protein=ProteinFactory(name=name), name="default")

    def names(query: str) -> list[str]:
        response = client.get(f"/api/proteins/spectra/?format=json&{query}")
        assert response.status_code == 200, response.content
        return sorted(p["name"] for p in response.json())

    assert names("name=spectraone") == ["SpectraOne"]
    assert names("name__icontains=spectrat") == ["SpectraThree", "SpectraTwo"]
    assert len(names("limit=2")) == 2
    assert len(names("")) == 3
    assert names("slug=spectratwo") == ["SpectraTwo"]
    assert client.get("/api/proteins/spectra/?format=json&protein=x").status_code == 400


@pytest.mark.django_db
def test_protein_list_api_search_alias(client):
    """`search=` is a name/alias search, as the conventional name for one."""
    ProteinFactory(name="GuessProtein")
    ProteinFactory(name="OtherProtein")
    response = client.get("/api/proteins/?format=json&search=guess")
    assert [p["name"] for p in response.json()] == ["GuessProtein"]


@pytest.mark.parametrize(
    ("param", "hint"),
    [
        ("pdb_id", "pdb"),
        ("pdb__icontains", "pdb__contains"),
        ("ex_max__gt", "ex_max__gte"),
        ("name__contains", "name__icontains"),
        ("default_state__ex_max__exact", "default_state__ex_max"),
        ("ex_maxx", "ex_max"),
        ("q", None),
        ("find", None),
    ],
)
@pytest.mark.django_db
def test_unknown_param_did_you_mean(client, param, hint):
    """A guessed param's 400 names the param it most likely meant, if any."""
    response = client.get(f"/api/proteins/?format=json&{param}=1")
    assert response.status_code == 400
    error = response.json()
    assert error["did_you_mean"] == ({param: hint} if hint else {})
    assert (f"Did you mean: {hint} (for {param})" in error["detail"]) == bool(hint)


@pytest.mark.django_db
def test_protein_list_api_ex_em_max_aliases(client):
    ProteinFactory(name="Green", default_state__ex_max=488, default_state__em_max=507)
    ProteinFactory(name="Red", default_state__ex_max=587, default_state__em_max=610)
    for query in ("ex_max=488", "em_max=507"):
        response = client.get(f"/api/proteins/?format=json&{query}")
        assert [p["name"] for p in response.json()] == ["Green"], query


@pytest.mark.django_db
def test_protein_spectra_api_query_count(client, django_assert_max_num_queries):
    """Spectra are prefetched, not queried per state."""
    for i in range(4):
        StateFactory(protein=ProteinFactory(name=f"Many{i}"), name="default")
    with django_assert_max_num_queries(8):
        response = client.get("/api/proteins/spectra/?format=json")
    assert len(response.json()) == 4
    assert all(p["spectra"] for p in response.json())


@pytest.mark.django_db
def test_protein_detail_by_other_identifiers(client):
    protein = ProteinFactory(name="Lookup Me", aliases=["LkM"], pdb=["1ABC", "2SHR"])
    other = ProteinFactory(name="Other", pdb=["2SHR"])
    for key in (protein.uuid.lower(), "Lookup%20Me", "lkm", "1abc"):
        response = client.get(f"/api/proteins/{key}/?format=json&fields=slug")
        assert response.json() == {"slug": protein.slug}, key

    # a PDB ID shared by two proteins: say which, rather than pick one
    response = client.get("/api/proteins/2SHR/?format=json")
    assert response.status_code == 404
    assert other.slug in response.json()["detail"] and "?pdb=2SHR" in response.json()["detail"]

    response = client.get("/api/proteins/no-such-thing/?format=json")
    assert response.status_code == 404
    assert "PDB ID" in response.json()["detail"]

    # an exact alias is found among many partial matches, whatever their order
    for i in range(30):
        ProteinFactory(name=f"AAA-LkM-{i}")
    assert client.get("/api/proteins/lkm/?format=json&fields=slug").json() == {
        "slug": protein.slug
    }

    # part of a name or alias is not a match (`cherry` is not PAmCherry1)
    for key in ("lookup", "km", "ookup%20m", "a%00b"):
        assert client.get(f"/api/proteins/{key}/?format=json").status_code == 404, key


@pytest.mark.django_db
def test_protein_list_page_and_page_size(client):
    for i in range(7):
        ProteinFactory(name=f"Page{i}")
    slugs = [p["slug"] for p in client.get("/api/proteins/?format=json").json()]

    def page(query):
        return [p["slug"] for p in client.get(f"/api/proteins/?format=json&{query}").json()]

    assert page("page=2&page_size=3") == slugs[3:6]
    assert page("page=3&page_size=3") == slugs[6:]
    assert page("page_size=2") == slugs[:2]
    assert page("page=1") == slugs  # (fewer than the default page size)
    assert page("page=x&page_size=y") == slugs[:100]  # (garbage is ignored, as for limit)


@pytest.mark.django_db
def test_protein_list_fields_and_include_fields(client):
    protein = ProteinFactory(name="Fields")
    state = protein.states.get()
    spectrum = state.spectra.get(subtype="ex")

    response = client.get("/api/proteins/?format=json&fields=name,states__ex_max")
    assert response.json() == [{"name": "Fields", "states": [{"ex_max": state.ex_max}]}]

    # spectra are on demand: left out of the full output, unless asked for by name
    url = "/api/proteins/fields/?format=json"
    assert "ex_spectrum" not in client.get(url).json()["states"][0]
    data = [list(p) for p in spectrum.data]
    response = client.get(f"{url}&include_fields=states__ex_spectrum")
    assert response.json()["states"][0]["ex_spectrum"] == data
    assert response.json()["states"][0]["ex_max"] == state.ex_max
    response = client.get(f"{url}&fields=slug,states__ex_spectrum")
    assert response.json() == {"slug": "fields", "states": [{"ex_spectrum": data}]}

    # the spectra endpoint builds its output from the states: fine without them
    response = client.get("/api/proteins/spectra/?format=json&fields=name")
    assert response.json() == [{"name": "Fields", "spectra": []}]
    # the table endpoint's serializer cannot choose fields: rejected, not ignored
    response = client.get("/api/proteins/table-data/?format=json&fields=name")
    assert response.status_code == 400
    assert "fields" in response.json()["detail"]


@pytest.mark.django_db
def test_protein_list_short_filter_names(client):
    ProteinFactory(name="Red", default_state__ex_max=590, default_state__em_max=610)
    ProteinFactory(name="Green", default_state__ex_max=488, default_state__em_max=510)

    def names(query):
        return [p["name"] for p in client.get(f"/api/proteins/?format=json&{query}").json()]

    assert names("ex_max__gte=500") == ["Red"]
    assert names("em_max__range=500,520") == ["Green"]
    assert names("ex_max__around=590") == ["Red"]
    assert names("default_state__ex_max__gte=500") == ["Red"]  # (still works)


@pytest.mark.django_db
def test_api_schema_and_docs_are_public(client):
    assert client.get("/api/schema/").status_code == 200
    assert client.get("/api/docs/").status_code == 200
    response = client.get("/api/openapi.json")
    assert response.status_code == 200
    schema = response.json()
    parameters = {p["name"] for p in schema["paths"]["/api/proteins/"]["get"]["parameters"]}
    assert {"fields", "include_fields", "page", "page_size", "limit", "offset"} <= parameters
    detail = schema["paths"]["/api/proteins/{slug}/"]["get"]["parameters"]
    assert "PDB ID" in next(p for p in detail if p["name"] == "slug")["description"]
