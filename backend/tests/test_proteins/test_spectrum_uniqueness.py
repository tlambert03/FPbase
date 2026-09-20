"""A fluorophore has at most one *approved* spectrum per subtype."""

from __future__ import annotations

import json

import pytest
from django.contrib.auth.models import Permission
from django.core.exceptions import ValidationError
from django.db import IntegrityError, transaction
from django.urls import reverse

from proteins.factories import (
    DyeStateFactory,
    FilterFactory,
    ProteinFactory,
    SpectrumFactory,
    StateFactory,
)
from proteins.forms.spectrum_v2 import SpectrumFormV2
from proteins.models import FluorState, Spectrum, State
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


def test_db_rejects_second_approved_spectrum_of_same_subtype():
    state = StateFactory()  # comes with approved ex, em, and 2p spectra
    dup = SpectrumFactory.build(owner_fluor=state, category="p", subtype="ex")
    with pytest.raises(ValidationError, match="already has an approved spectrum"):
        dup.full_clean()
    with pytest.raises(IntegrityError), transaction.atomic():
        dup.save()
    assert state.ex_spectrum is not None


@pytest.mark.parametrize("status", ["pending", "rejected"])
def test_unapproved_spectra_may_coexist_with_approved(status):
    state = StateFactory()
    for _ in range(2):
        dup = SpectrumFactory.build(owner_fluor=state, category="p", subtype="ex", status=status)
        dup.full_clean()
        dup.save()
    assert state.ex_spectrum is not None


def test_other_subtype_still_allowed():
    state = StateFactory()
    SpectrumFactory(owner_fluor=state, category="p", subtype="ab")
    assert state.spectra.count() == 4


def _form(category: str, owner: str, subtype: str, owner_slug: str | None = None):
    spec = {
        "data": [[400, 0.1], [500, 1.0], [600, 0.5]],
        "category": category,
        "owner": owner,
        "owner_slug": owner_slug,
        "subtype": subtype,
        "scale_factor": None,
        "ph": None,
        "solvent": None,
        "peak_wave": 500,
        "column_name": "A",
    }
    data = {"spectra_json": json.dumps([spec]), "source": "test", "confirmation": True}
    return SpectrumFormV2(data)


def test_v2_form_rejects_existing_approved_protein_spectrum():
    state = StateFactory()
    form = _form("p", state.protein.name, "ex", owner_slug=state.protein.slug)
    assert not form.is_valid()
    assert "already has a spectrum" in str(form.errors["spectra_json"])


def test_v2_form_allows_resubmission_while_pending():
    state = StateFactory(ex_spectrum=None)
    assert not state.spectra.filter(subtype="ex").exists()
    for _ in range(2):  # e.g. a submitter correcting their own pending submission
        form = _form("p", state.protein.name, "ex", owner_slug=state.protein.slug)
        assert form.is_valid(), form.errors
        form.save()
    pending = Spectrum.objects.all_objects().filter(owner_fluor=state, subtype="ex")
    assert [s.status for s in pending] == ["pending", "pending"]


def test_v2_form_rejects_existing_dye_and_filter_spectra():
    dye_state = DyeStateFactory(name=FluorState.DEFAULT_NAME)  # the state the form submits to
    assert not _form("d", dye_state.dye.name, "em").is_valid()

    spectrum = FilterFactory().spectrum
    assert not _form("f", spectrum.owner_filter.name, spectrum.subtype).is_valid()


def test_v2_form_accepts_new_subtype():
    state = StateFactory()
    form = _form("p", state.protein.name, "ab", owner_slug=state.protein.slug)
    assert form.is_valid(), form.errors
    (created,) = form.save()
    assert created.owner_fluor_id == state.id


@pytest.fixture
def moderator_client(client):
    moderator = UserFactory()
    perms = Permission.objects.filter(codename__in=["change_spectrum", "delete_spectrum"])
    moderator.user_permissions.set(perms)
    client.force_login(moderator)
    return client


def _accept(client, *spectra: Spectrum):
    url = reverse("proteins:pending_spectrum_action")
    return client.post(url, {"spectrum_ids[]": [s.id for s in spectra], "action": "accept"})


def test_cannot_accept_pending_spectrum_when_approved_exists(moderator_client):
    state = StateFactory()
    pending = SpectrumFactory(owner_fluor=state, category="p", subtype="ex", status="pending")

    response = _accept(moderator_client, pending)
    assert response.status_code == 409
    assert "Delete the approved spectrum first" in response.json()["error"]
    pending.refresh_from_db()
    assert pending.status == "pending"

    state.ex_spectrum.delete()
    assert _accept(moderator_client, pending).json()["success"]
    assert State.objects.get(id=state.id).ex_spectrum == pending


def test_cannot_accept_two_pending_duplicates_at_once(moderator_client):
    state = StateFactory(ex_spectrum=None)
    dups = [
        SpectrumFactory(owner_fluor=state, category="p", subtype="ex", status="pending")
        for _ in range(2)
    ]
    assert _accept(moderator_client, *dups).status_code == 409
    assert _accept(moderator_client, dups[0]).json()["success"]


def test_v2_form_submits_to_default_state():
    protein = ProteinFactory(default_state=None)
    StateFactory(protein=protein, name="first")
    default = StateFactory(protein=protein, name="second")
    protein.default_state = default
    protein.save()

    form = _form("p", protein.name, "ab", owner_slug=protein.slug)
    assert form.is_valid(), form.errors
    (created,) = form.save()
    assert created.owner_fluor_id == default.id
