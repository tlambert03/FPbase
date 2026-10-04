from __future__ import annotations

from importlib import import_module

import pytest
from django.apps import apps

from proteins.factories import FilterFactory
from proteins.models import Filter

pytestmark = pytest.mark.django_db

migration = import_module("proteins.migrations.0067_semrock_product_links")
SEARCH = "https://www.idex-hs.com/store/search-results/1/?searchCriteria="


def test_semrock_filter_links_to_a_store_search():
    semrock = FilterFactory(manufacturer="Semrock", part="FF01-520/35")
    assert semrock.url == SEARCH + "FF01-520/35"


def test_migration_rewrites_dead_semrock_links():
    semrock = FilterFactory(manufacturer="Semrock", part="FF01-520/35")
    chroma = FilterFactory(manufacturer="Chroma", part="ET525/50m")
    # as stored for filters saved before 2019: a dead link and a slug without the dash
    Filter.objects.filter(pk=semrock.pk).update(
        url=migration.OLD + "FF01-520/35", slug="semrock-ff01-52035"
    )

    migration.forward(apps, None)
    semrock.refresh_from_db()
    assert semrock.url == SEARCH + "FF01-520/35"
    assert semrock.slug == "semrock-ff01-52035"
    assert Filter.objects.get(pk=chroma.pk).url == chroma.url

    migration.backward(apps, None)
    semrock.refresh_from_db()
    assert semrock.url == migration.OLD + "FF01-520/35"
