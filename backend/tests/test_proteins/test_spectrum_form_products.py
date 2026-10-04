"""Filters, cameras and lights submitted through the spectrum form get a manufacturer."""

from __future__ import annotations

import json

import pytest

from proteins.factories import FilterFactory
from proteins.forms.spectrum_v2 import SpectrumFormV2
from proteins.models import Filter

pytestmark = pytest.mark.django_db


def _submit(owner: str, category: str = "f", subtype: str = "bp", **fields: str | None):
    spec = {
        "data": [[400, 0.1], [500, 0.9], [600, 0.1]],
        "category": category,
        "owner": owner,
        "owner_slug": None,
        "subtype": subtype,
        "scale_factor": None,
        "ph": None,
        "solvent": None,
        "peak_wave": None,
        "column_name": "A",
        **fields,
    }
    form = SpectrumFormV2(
        {"spectra_json": json.dumps([spec]), "source": "s", "confirmation": True}
    )
    assert form.is_valid(), form.errors
    (spectrum,) = form.save()
    return spectrum.owner


def test_manufacturer_and_part_are_saved():
    filt = _submit("My filter", manufacturer="Semrock", part="FF01-520/35")
    assert (filt.manufacturer, filt.part) == ("Semrock", "FF01-520/35")
    assert filt.url == "https://www.idex-hs.com/store/search-results/1/?searchCriteria=FF01-520/35"


def test_manufacturer_takes_the_existing_spelling():
    FilterFactory(manufacturer="Thorlabs")
    assert _submit("A Thorlabs filter", manufacturer="ThorLabs").manufacturer == "Thorlabs"


def test_manufacturer_is_inferred_from_the_name():
    FilterFactory(manufacturer="Chroma")
    assert _submit("Chroma ET525/50m").manufacturer == "Chroma"
    assert _submit("chroma-ET470/40x").manufacturer == "Chroma"
    assert _submit("Chromatic aberration filter").manufacturer == ""
    assert _submit("My custom dichroic").manufacturer == ""


def test_new_manufacturer_is_kept_as_typed():
    assert _submit("X", manufacturer="  Acme  Optics ").manufacturer == "Acme Optics"


def test_owner_name_whitespace_is_collapsed():
    assert _submit("  Semrock   FF01-520/35 ").name == "Semrock FF01-520/35"


def test_light_gets_a_manufacturer_too():
    light = _submit("Some LED", category="l", subtype="pd", manufacturer="CoolLED")
    assert light.manufacturer == "CoolLED"


def test_fluorophores_ignore_product_fields():
    dye_state = _submit("Some dye", category="d", subtype="em", manufacturer="Biotium")
    assert dye_state.dye.manufacturer == ""
    assert not Filter.objects.exists()
