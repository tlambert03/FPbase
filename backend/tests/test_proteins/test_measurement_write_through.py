"""Direct edits of a state's fluorescence values are recorded as measurements.

Measurements are the source of truth and the state's values are a summary of them, so an
edit that only touched the state would be lost by the next `rebuild_attributes()`.
"""

from __future__ import annotations

from typing import cast

import pytest
import reversion
from django.core.management import call_command
from django.urls import reverse
from reversion.models import Version

from proteins.models import Dye, DyeState, Protein, State
from proteins.models import FluorescenceMeasurement as FM
from references.models import Reference
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


def _ref(doi: str) -> Reference:
    ref = Reference(doi=doi, year=2020)
    ref.save(skipdoi=True)
    return ref


@pytest.fixture
def protein() -> Protein:
    return Protein.objects.create(name="TestProtein", primary_reference=_ref("10.1234/primary"))


@pytest.fixture
def state(protein: Protein) -> State:
    return State.objects.create(protein=protein, name="default", ex_max=488, qy=0.6)


def test_new_state_values_become_a_primary_reference_measurement(state: State):
    (m,) = state.measurements.all()
    assert (m.ex_max, m.qy, m.em_max) == (488, 0.6, None)
    assert m.reference == state.protein.primary_reference

    state.rebuild_attributes()
    state.refresh_from_db()
    assert (state.ex_max, state.qy) == (488, 0.6)
    assert state.source_map["ex_max"] == m.id


def test_state_without_values_gets_no_measurement(protein: Protein):
    state = State.objects.create(protein=protein, name="empty")
    assert not state.measurements.exists()


def test_edit_updates_the_primary_measurement(state: State):
    state = State.objects.get(id=state.id)
    state.ex_max = 490
    state.save()

    (m,) = state.measurements.all()
    assert (m.ex_max, m.qy) == (490, 0.6)
    state.rebuild_attributes()
    assert State.objects.get(id=state.id).ex_max == 490


def test_clearing_a_value_clears_the_measurement(state: State):
    state = State.objects.get(id=state.id)
    state.qy = None
    state.save()
    state.rebuild_attributes()
    assert State.objects.get(id=state.id).qy is None


def test_unchanged_save_writes_nothing(state: State):
    # a value that came from another paper must not be copied to the primary measurement
    FM.objects.create(state=state, reference=_ref("10.1234/other"), pka=5.5)
    state = State.objects.get(id=state.id)
    assert state.pka == 5.5
    modified = list(FM.objects.order_by("id").values_list("modified", flat=True))

    state.name = "renamed"
    state.save()
    assert list(FM.objects.order_by("id").values_list("modified", flat=True)) == modified
    assert state.measurements.get(reference=state.protein.primary_reference).pka is None


def test_edit_outranks_a_value_from_another_reference(state: State):
    FM.objects.create(state=state, reference=_ref("10.1234/other"), pka=5.5)
    state = State.objects.get(id=state.id)
    state.pka = 6.0
    state.save()

    state = State.objects.get(id=state.id)
    assert state.pka == 6.0
    assert state.measurements.get(reference=state.protein.primary_reference).pka == 6.0
    assert state.measurements.get(reference__doi="10.1234/other").pka == 5.5


def test_edit_clears_pin_on_that_field(state: State):
    other = FM.objects.create(state=state, reference=_ref("10.1234/other"), qy=0.1, pka=5.5)
    state = State.objects.get(id=state.id)
    state.pinned_source_map = {"qy": other.id, "pka": other.id}
    state.save()
    state.rebuild_attributes()
    assert state.qy == 0.1

    state.qy = 0.7
    state.save()
    state = State.objects.get(id=state.id)
    assert state.qy == 0.7
    assert state.pinned_source_map == {"pka": other.id}


def test_stale_state_is_not_written_back(protein: Protein):
    # rebuild_cache=False leaves the state stale; saving it must not clobber the evidence
    state = State.objects.create(protein=protein, name="default")
    FM(state=state, reference=protein.primary_reference, ex_max=488).save(rebuild_cache=False)
    state = State.objects.get(id=state.id)
    assert state.ex_max is None

    state.name = "renamed"
    state.save()
    assert state.measurements.get().ex_max == 488


def test_dye_edits_go_to_an_unattributed_measurement():
    dye = Dye.objects.create(name="TestDye")
    dye_state = DyeState.objects.create(dye=dye, name="default", ex_max=550)
    (m,) = dye_state.measurements.all()
    assert m.reference is None

    # with no primary reference, the unattributed measurement wins over newer ones
    FM.objects.create(state=dye_state, reference=_ref("10.1234/other"), ex_max=560)
    assert DyeState.objects.get(id=dye_state.id).ex_max == 550


def test_primary_reference_change_keeps_values_and_attributes_new_edits(state: State):
    protein = state.protein
    old_ref, new_ref = protein.primary_reference, _ref("10.1234/new")
    protein.primary_reference = new_ref
    protein.save()

    state = protein.states.get()
    state.ex_max = 495
    state.save()

    state.rebuild_attributes()
    state = State.objects.get(id=state.id)
    assert (state.ex_max, state.qy) == (495, 0.6)
    assert state.measurements.get(reference=new_ref).ex_max == 495
    assert state.measurements.get(reference=old_ref).ex_max == 488


def test_protein_form_submission_and_edit(client):
    user = UserFactory()
    client.force_login(user)
    data = {
        "name": "FormProtein",
        "reference_doi": "10.1038/nmeth.2413",
        "states-0-name": "default",
        "states-0-ex_max": 488,
        "states-0-em_max": 525,
        "confirmation": True,
        "lineage-TOTAL_FORMS": 1,
        "lineage-INITIAL_FORMS": 0,
        "lineage-MIN_NUM_FORMS": 0,
        "lineage-MAX_NUM_FORMS": 1,
        "states-TOTAL_FORMS": 1,
        "states-INITIAL_FORMS": 0,
        "states-MIN_NUM_FORMS": 0,
        "states-MAX_NUM_FORMS": 1000,
    }
    assert client.post(reverse("proteins:submit"), data).status_code == 302

    protein = Protein.objects.get(name="FormProtein")
    state = cast("State", protein.default_state)
    (m,) = state.measurements.all()
    assert (m.ex_max, m.em_max) == (488, 525)
    assert m.reference == protein.primary_reference
    assert m.created_by == user

    data |= {
        "states-INITIAL_FORMS": 1,
        "states-0-fluorstate_ptr": state.id,
        "states-0-ex_max": 490,
    }
    url = reverse("proteins:update", args=[protein.slug])
    assert client.post(url, data).status_code == 302

    state.rebuild_attributes()
    state = State.objects.get(id=state.id)
    assert (state.ex_max, state.em_max) == (490, 525)
    assert state.measurements.count() == 1


def test_revision_rollback_survives_rebuild(client, protein: Protein):
    with reversion.create_revision():
        protein.save()
        state = State.objects.create(protein=protein, name="default", ex_max=488)
    with reversion.create_revision():
        state.ex_max = 500
        state.save()

    # (reversion restores rows with save_base(), bypassing State.save)
    client.force_login(UserFactory(is_staff=True))
    revision = Version.objects.get_for_object(state).last().revision
    url = reverse("proteins:admin_revert_revision", args=[revision.id])
    assert client.post(url).status_code == 200

    state = State.objects.get(id=state.id)
    assert state.ex_max == 488
    state.rebuild_attributes()
    assert State.objects.get(id=state.id).ex_max == 488


def test_sync_measurements_command(protein: Protein, capsys):
    # states written before write-through existed: values with no (or stale) measurements
    unmeasured = State.objects.create(protein=protein, name="unmeasured")
    edited = State.objects.create(protein=protein, name="edited", ex_max=488)
    State.objects.filter(id=unmeasured.id).update(ex_max=400, qy=0.5)
    State.objects.filter(id=edited.id).update(ex_max=490)

    call_command("sync_measurements", "--dry-run")
    out = capsys.readouterr().out
    assert "unmeasured" in out and "edited" in out
    assert not unmeasured.measurements.exists()

    call_command("sync_measurements")
    for state, expected in ((unmeasured, 400), (edited, 490)):
        state.rebuild_attributes()
        assert State.objects.get(id=state.id).ex_max == expected
    assert edited.measurements.count() == 1

    call_command("sync_measurements")
    assert "0 state(s)" in capsys.readouterr().out
