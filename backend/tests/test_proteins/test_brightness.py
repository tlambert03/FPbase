"""Brightness is derived from extinction coefficient and quantum yield."""

from __future__ import annotations

import pytest

from proteins.models import FluorescenceMeasurement as FM
from proteins.models import Protein, State

pytestmark = pytest.mark.django_db


@pytest.fixture
def state() -> State:
    protein = Protein.objects.create(name="TestProtein")
    return State.objects.create(protein=protein, name="default", ext_coeff=100_000, qy=0.5)


def test_brightness_is_computed(state: State):
    assert state.brightness == 50


@pytest.mark.parametrize(
    ("field", "value", "expected"), [("qy", 0, 0), ("qy", None, None), ("ext_coeff", None, None)]
)
def test_brightness_follows_its_inputs(state: State, field: str, value, expected):
    setattr(state, field, value)
    state.save()
    state.refresh_from_db()
    assert state.brightness == expected


def test_measurement_brightness_follows_its_inputs(state: State):
    m = FM.objects.create(state=state, ext_coeff=100_000, qy=0.5)
    assert m.brightness == 50
    m.qy = 0
    m.save()
    assert m.brightness == 0
