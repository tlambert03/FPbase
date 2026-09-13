from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model
from django.core.cache import cache

from favit.models import Favorite
from proteins import search_index
from proteins.extrest.ga import ga_spectra_views
from proteins.factories import DyeFactory, ProteinFactory
from proteins.models import Protein, Spectrum

URL = "/api/search-index/"


@pytest.fixture(autouse=True)
def _clear_cache(monkeypatch: pytest.MonkeyPatch):
    cache.clear()
    monkeypatch.setattr(search_index, "cached_ga_popular", lambda: {"year": []})
    monkeypatch.setattr(search_index, "cached_ga_spectra_views", dict)
    yield
    cache.clear()


def _by_slug(data: dict) -> dict[str, dict]:
    return {p["slug"]: p for p in data["proteins"]}


@pytest.mark.django_db
def test_search_index_payload(client) -> None:
    visible = ProteinFactory(name="mTestVisible", aliases=["TV1"], genbank="AB123")
    ProteinFactory(name="mTestHidden", status=Protein.STATUS.hidden)

    response = client.get(URL)
    assert response.status_code == 200
    assert response["Content-Type"] == "application/json"
    assert "public" in response["Cache-Control"]
    data = response.json()

    assert data["version"] == search_index.SEARCH_INDEX_VERSION
    prots = _by_slug(data)
    assert "mtesthidden" not in prots
    rec = prots[visible.slug]
    assert rec["name"] == "mTestVisible"
    assert rec["url"] == visible.get_absolute_url()
    assert rec["aliases"] == ["TV1"]
    assert rec["genbank"] == "AB123"
    assert rec["spectra"] is True  # StateFactory adds ex/em spectra
    assert rec["color"] == visible.color

    ref = next(
        r for r in data["references"] if r["url"] == visible.primary_reference.get_absolute_url()
    )
    assert "mTestVisible" in ref["primary"]
    assert data["organisms"][0]["url"].startswith("/organism/")


@pytest.mark.django_db
def test_search_index_unescapes_html_entities(client) -> None:
    ProteinFactory(name="mKO&kappa;", slug="mkokappa")
    assert _by_slug(client.get(URL).json())["mkokappa"]["name"] == "mKOκ"


@pytest.mark.django_db
def test_search_index_popularity(client, monkeypatch: pytest.MonkeyPatch) -> None:
    viewed, faved, obscure = (ProteinFactory(name=n) for n in ("mViewed", "mFaved", "mObscure"))
    views = {"year": [(viewed.slug, viewed.name, 90.0), (obscure.slug, obscure.name, 0.001)]}
    monkeypatch.setattr(search_index, "cached_ga_popular", lambda: views)
    User = get_user_model()
    for i in range(3):
        user = User.objects.create_user(username=f"user{i}", password="pw")
        Favorite.objects.create(user, faved.id, "proteins.Protein")

    prots = _by_slug(client.get(URL).json())
    assert prots[viewed.slug]["p"] == 1
    assert 0 < prots[faved.slug]["p"] < 1
    assert "p" not in prots[obscure.slug] or prots[obscure.slug]["p"] < prots[faved.slug]["p"]


@pytest.mark.django_db
def test_search_index_ga_failure_falls_back(client, monkeypatch: pytest.MonkeyPatch) -> None:
    def boom():
        raise ValueError("no credentials")

    monkeypatch.setattr(search_index, "cached_ga_popular", boom)
    monkeypatch.setattr(search_index, "cached_ga_spectra_views", boom)
    ProteinFactory()
    DyeFactory()
    assert client.get(URL).status_code == 200


@pytest.mark.django_db
def test_search_index_etag(client) -> None:
    protein = ProteinFactory()
    response = client.get(URL)
    etag = response["ETag"]
    assert etag.startswith('W/"')
    assert client.get(URL, headers={"If-None-Match": etag}).status_code == 304

    protein.name = "mRenamed"
    protein.save()
    response = client.get(URL, headers={"If-None-Match": etag})
    assert response.status_code == 200
    assert response["ETag"] != etag
    assert "mRenamed" in {p["name"] for p in response.json()["proteins"]}


@pytest.mark.django_db
def test_search_index_query_count(django_assert_max_num_queries) -> None:
    for _ in range(10):
        ProteinFactory()
    with django_assert_max_num_queries(12):
        search_index.build_search_index()


@pytest.mark.django_db
def test_search_index_dyes(client, monkeypatch: pytest.MonkeyPatch) -> None:
    popular, obscure = DyeFactory(name="Alexa Fluor 488"), DyeFactory(name="Obscure Dye")
    ex_em = set(
        Spectrum.objects.filter(
            owner_fluor__dyestate__dye=popular, subtype__in=("ex", "em")
        ).values_list("id", flat=True)
    )
    monkeypatch.setattr(search_index, "cached_ga_spectra_views", lambda: {min(ex_em): 50})

    dyes = {d["name"]: d for d in client.get(URL).json()["dyes"]}
    rec = dyes[popular.name]
    path, ids = rec["url"].split("?s=")
    assert path == "/spectra/"
    assert {int(i) for i in ids.split(",")} == ex_em  # 2p spectra are left out
    assert rec["p"] == 1
    assert dyes[obscure.name].get("p", 0) < 1


class _Value:
    def __init__(self, value: str) -> None:
        self.value = value


class _Row:
    def __init__(self, path: str, views: int) -> None:
        self.dimension_values = [_Value(path)]
        self.metric_values = [_Value(str(views))]


def test_ga_spectra_views_parses_viewer_urls() -> None:
    rows = [
        _Row("/spectra/?s=17,18,$cl0_488&showY=0", 10),
        _Row("/spectra/?xMin=400&s=18,18", 5),
        _Row("/spectra/?showY=1", 99),
    ]
    client = type("Client", (), {"run_report": lambda self, req: type("R", (), {"rows": rows})})()
    assert ga_spectra_views(client) == {17: 10, 18: 15}
