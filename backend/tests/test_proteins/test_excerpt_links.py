"""Protein names in excerpts are linked to their pages."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from django.core.cache import cache

from proteins.factories import ProteinFactory
from proteins.util.helpers import link_excerpts

pytestmark = pytest.mark.django_db


@pytest.fixture(autouse=True)
def _clear_cache():
    cache.clear()
    yield
    cache.clear()


def _linked(content: str, **kwargs) -> str:
    excerpts = link_excerpts([SimpleNamespace(content=content)], **kwargs)
    assert excerpts
    return excerpts[0].content


def _link(protein) -> str:
    return f'<a href="{protein.get_absolute_url()}" class="text-info">{protein.name}</a>'


def test_other_proteins_are_linked():
    protein = ProteinFactory(name="mFoo")
    assert _linked("We compared mFoo to others") == f"We compared {_link(protein)} to others"


def test_own_name_and_aliases_are_emphasized():
    ProteinFactory(name="mFoo", aliases=["Foo2"])
    out = _linked("Both mFoo and Foo2 here", obj_name="mFoo", aliases=["Foo2"])
    assert out == "Both <strong>mFoo</strong> and <strong>Foo2</strong> here"


def test_names_are_matched_literally():
    star, paren = ProteinFactory(name="Padron*"), ProteinFactory(name="GFP (S65T)")
    out = _linked("Using Padron and Padron* with GFP (S65T) here")
    assert out == f"Using Padron and {_link(star)} with {_link(paren)} here"


def test_content_is_shown_as_written():
    protein = ProteinFactory(name="mFoo")
    out = _linked("Cycles <10 for mFoo & <b>others</b>")
    assert out == f"Cycles &lt;10 for {_link(protein)} &amp; &lt;b&gt;others&lt;/b&gt;"


def test_names_are_shown_as_written():
    protein = ProteinFactory(name="mFoo<i>&")
    link = f'<a href="{protein.get_absolute_url()}" class="text-info">mFoo&lt;i&gt;&amp;</a>'
    assert _linked("Using mFoo<i>& here") == f"Using {link} here"
