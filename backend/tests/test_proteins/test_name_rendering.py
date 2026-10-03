"""Protein names are displayed as written, whatever characters they contain."""

from __future__ import annotations

import json
import re

import pytest
from django.core.cache import cache
from django.urls import reverse
from django.utils.html import escape

from proteins.factories import ProteinFactory
from proteins.models import Protein

pytestmark = pytest.mark.django_db

NAME = '<b>FP</b> & "co"'


@pytest.fixture(autouse=True)
def _clear_cache():
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def protein() -> Protein:
    return ProteinFactory(name=NAME)


def _assert_name_escaped(content: str) -> None:
    assert escape(NAME) in content
    assert NAME not in content
    assert "<b>FP</b>" not in content


def test_protein_detail(client, protein: Protein):
    _assert_name_escaped(client.get(protein.get_absolute_url()).content.decode())


def test_protein_detail_structured_data(client, protein: Protein):
    content = client.get(protein.get_absolute_url()).content.decode()
    blocks = re.findall(r'<script type="application/ld\+json">(.*?)</script>', content, re.DOTALL)
    dataset = next(d for d in map(json.loads, blocks) if d["@type"] == "Dataset")
    assert dataset["name"] == NAME
    assert dataset["description"].startswith(NAME)


def test_structured_data_leaves_out_missing_values(client):
    protein = ProteinFactory(name="chromo")
    state = protein.default_state
    state.em_max = state.qy = None
    state.save()
    content = client.get(protein.get_absolute_url()).content.decode()
    blocks = re.findall(r'<script type="application/ld\+json">(.*?)</script>', content, re.DOTALL)
    dataset = next(d for d in map(json.loads, blocks) if d["@type"] == "Dataset")
    measured = {v["name"]: v["value"] for v in dataset["variableMeasured"]}
    assert measured["Excitation Maximum"] == protein.default_state.ex_max
    assert "Emission Maximum" not in measured
    assert "Quantum Yield" not in measured


def test_activity(client, protein: Protein):
    _assert_name_escaped(client.get(reverse("proteins:activity")).content.decode())


def _fake_muscle_html(queryset, output: str) -> tuple[str, str]:
    """Stand-in for muscle's html alignment: 3 header lines, then a row per sequence."""
    rows = [f"{p.uuid:16}<span style=background-color:#FFEEE0>{p.seq}</span>" for p in queryset]
    return "\n".join(["", "", "", *rows]), ""


@pytest.mark.parametrize("n_others", [1, 2])
def test_compare(client, protein: Protein, n_others: int, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(type(Protein.objects.all()), "to_tree", _fake_muscle_html)
    slugs = [protein.slug, *(p.slug for p in ProteinFactory.create_batch(n_others))]
    url = reverse("proteins:compare", args=(",".join(slugs),))
    content = client.get(url).content.decode()
    _assert_name_escaped(content)
    assert ("<span style=background-color:#FFFFFF>" in content) == (n_others > 1)
