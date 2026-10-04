from __future__ import annotations

import pytest

from proteins.factories import DyeFactory, FilterFactory, SpectrumFactory, StateFactory
from proteins.models import Spectrum
from proteins.models.spectrum import get_spectra_list

pytestmark = pytest.mark.django_db

THERMO = "https://www.thermofisher.com/order/catalog/product/A10169"


def _owner_url(spectrum: Spectrum) -> str | None:
    (item,) = (s for s in get_spectra_list() if s["id"] == spectrum.id)
    return item["owner"]["url"]


def test_protein_spectrum_links_to_the_protein():
    state = StateFactory()
    spectrum = state.spectra.first()
    assert _owner_url(spectrum) == state.protein.slug


def test_dye_spectrum_links_to_the_vendor():
    linked = DyeFactory(url=THERMO).states.get()
    unlinked = DyeFactory().states.get()
    spectra = [
        SpectrumFactory(owner_fluor=state, category=Spectrum.DYE, subtype=Spectrum.ABS)
        for state in (linked, unlinked)
    ]
    assert _owner_url(spectra[0]) == THERMO
    assert not _owner_url(spectra[1])


def test_filter_spectrum_links_to_the_vendor():
    filt = FilterFactory(manufacturer="Chroma", part="ET525/50m")
    assert _owner_url(filt.spectrum) == "https://www.chroma.com/products/parts/ET525-50m"
