"""Bulk optical-config entry on the microscope form looks filters up in the database."""

from __future__ import annotations

import pytest

from proteins.factories import FilterFactory
from proteins.forms.microscope import MicroscopeForm
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db

CONFIG = "Widefield Green, ET470/40x, T495lpxr, ET525/50m"


def _form(optical_configs: str) -> MicroscopeForm:
    data = {"name": "scope", "optical_configs": optical_configs}
    return MicroscopeForm(data, user=UserFactory())


def test_bulk_config_finds_filters_by_name():
    names = ("Chroma ET470/40x", "Chroma T495lpxr", "Chroma ET525/50m")
    filters = {FilterFactory(name=name) for name in names}
    form = _form(CONFIG)
    assert form.is_valid(), form.errors
    (config,) = form.save().optical_configs.all()
    assert set(config.filters.all()) == filters


def test_bulk_config_with_an_unknown_filter_is_an_error():
    form = _form(CONFIG)
    assert not form.is_valid()
    assert "Filter not found in database: ET470/40x" in form.errors["optical_configs"]
