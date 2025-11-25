"""Tests for Fluorophore.rebuild_attributes() compositing logic."""

from datetime import date

import pytest

from proteins.models import Dye, DyeState, Protein, State
from proteins.models import FluorescenceMeasurement as FM
from references.models import Reference


@pytest.fixture
def protein(db) -> Protein:
    """Create a simple protein."""
    return Protein.objects.create(name="TestProtein")


@pytest.fixture
def protein_with_primary_ref(db) -> Protein:
    """Create a protein with a primary reference."""
    ref = Reference(doi="10.1234/primary", year=2020)
    ref.save(skipdoi=True)
    return Protein.objects.create(name="TestProteinWithRef", primary_reference=ref)


@pytest.fixture
def dye(db) -> Dye:
    """Create a simple dye."""
    return Dye.objects.create(name="TestDye")


@pytest.fixture
def dye_with_primary_ref(db) -> Dye:
    """Create a dye with a primary reference."""
    ref = Reference(doi="10.1234/dye-primary", year=2020)
    ref.save(skipdoi=True)
    return Dye.objects.create(name="TestDyeWithRef", primary_reference=ref)


@pytest.fixture
def state(protein: Protein) -> State:
    """Create a state for the protein."""
    return State.objects.create(protein=protein, name="default")


@pytest.fixture
def state_with_ref(protein_with_primary_ref: Protein) -> State:
    """Create a state for a protein with primary reference."""
    return State.objects.create(protein=protein_with_primary_ref, name="default")


@pytest.fixture
def dyestate(dye: Dye) -> DyeState:
    """Create a dyestate for the dye."""
    return DyeState.objects.create(dye=dye, name="default")


@pytest.fixture
def dyestate_with_ref(dye_with_primary_ref: Dye) -> DyeState:
    """Create a dyestate for a dye with primary reference."""
    return DyeState.objects.create(dye=dye_with_primary_ref, name="default")


class TestRebuildAttributesBasicWaterfall:
    """Test basic waterfall logic - first non-null value wins."""

    def test_single_measurement_sets_values(self, state: State):
        """A single measurement should populate the fluorophore."""
        m = FM(fluorophore=state, ex_max=488, em_max=509, qy=0.67)
        m.save(rebuild_cache=False)
        state.rebuild_attributes()
        state.refresh_from_db()

        assert state.ex_max == 488
        assert state.em_max == 509
        assert state.qy == 0.67

    def test_multiple_measurements_first_value_wins(self, state: State):
        """When multiple measurements exist, first non-null value wins."""
        # Measurement with partial data (no qy)
        m1 = FM(fluorophore=state, ex_max=488, em_max=509, qy=None, date_measured=date(2020, 1, 1))
        m1.save(rebuild_cache=False)
        # Measurement with different values
        m2 = FM(fluorophore=state, ex_max=500, em_max=520, qy=0.5, date_measured=date(2019, 1, 1))
        m2.save(rebuild_cache=False)

        state.rebuild_attributes()
        state.refresh_from_db()

        # More recent measurement wins for ex_max and em_max
        assert state.ex_max == 488
        assert state.em_max == 509
        # qy comes from second measurement (first had null)
        assert state.qy == 0.5

    def test_source_map_tracks_measurement_ids(self, state: State):
        """source_map should track which measurement each field came from."""
        m1 = FM(fluorophore=state, ex_max=488, qy=None, date_measured=date(2020, 1, 1))
        m1.save(rebuild_cache=False)
        m2 = FM(fluorophore=state, ex_max=None, qy=0.5, date_measured=date(2019, 1, 1))
        m2.save(rebuild_cache=False)

        state.rebuild_attributes()
        state.refresh_from_db()

        assert state.source_map["ex_max"] == m1.id
        assert state.source_map["qy"] == m2.id


class TestRebuildAttributesTrustedPriority:
    """Test that is_trusted=True measurements take highest priority."""

    def test_trusted_measurement_overrides_recent(self, state: State):
        """A trusted measurement should override a more recent one."""
        # Recent but not trusted
        m1 = FM(fluorophore=state, ex_max=500, em_max=520, is_trusted=False, date_measured=date(2023, 1, 1))
        m1.save(rebuild_cache=False)
        # Older but trusted
        m2 = FM(fluorophore=state, ex_max=488, em_max=509, is_trusted=True, date_measured=date(2010, 1, 1))
        m2.save(rebuild_cache=False)

        state.rebuild_attributes()
        state.refresh_from_db()

        # Trusted measurement wins
        assert state.ex_max == 488
        assert state.em_max == 509


class TestRebuildAttributesPrimaryReferencePriority:
    """Test that primary reference measurements take priority over non-primary."""

    def test_primary_ref_measurement_overrides_recent(self, state_with_ref: State):
        """Measurement from primary reference should override more recent non-primary."""
        primary_ref = state_with_ref.protein.primary_reference
        other_ref = Reference(doi="10.1234/other", year=2022)
        other_ref.save(skipdoi=True)

        # Recent but not from primary reference
        m1 = FM(
            fluorophore=state_with_ref,
            reference=other_ref,
            ex_max=500,
            em_max=520,
            date_measured=date(2023, 1, 1),
        )
        m1.save(rebuild_cache=False)
        # Older but from primary reference
        m2 = FM(
            fluorophore=state_with_ref,
            reference=primary_ref,
            ex_max=488,
            em_max=509,
            date_measured=date(2010, 1, 1),
        )
        m2.save(rebuild_cache=False)

        state_with_ref.rebuild_attributes()
        state_with_ref.refresh_from_db()

        # Primary reference measurement wins
        assert state_with_ref.ex_max == 488
        assert state_with_ref.em_max == 509

    def test_primary_ref_on_dye(self, dyestate_with_ref: DyeState):
        """Test that primary reference works for dyes too."""
        primary_ref = dyestate_with_ref.dye.primary_reference
        other_ref = Reference(doi="10.1234/dye-other", year=2022)
        other_ref.save(skipdoi=True)

        # Not from primary reference
        m1 = FM(
            fluorophore=dyestate_with_ref,
            reference=other_ref,
            ex_max=600,
            em_max=650,
            date_measured=date(2023, 1, 1),
        )
        m1.save(rebuild_cache=False)
        # From primary reference
        m2 = FM(
            fluorophore=dyestate_with_ref,
            reference=primary_ref,
            ex_max=550,
            em_max=580,
            date_measured=date(2010, 1, 1),
        )
        m2.save(rebuild_cache=False)

        dyestate_with_ref.rebuild_attributes()
        dyestate_with_ref.refresh_from_db()

        # Primary reference measurement wins
        assert dyestate_with_ref.ex_max == 550
        assert dyestate_with_ref.em_max == 580


class TestRebuildAttributesPriorityOrder:
    """Test the complete priority order: trusted > primary_ref > date."""

    def test_trusted_beats_primary_ref(self, state_with_ref: State):
        """is_trusted should override even primary reference."""
        primary_ref = state_with_ref.protein.primary_reference

        # From primary reference but not trusted
        m1 = FM(
            fluorophore=state_with_ref,
            reference=primary_ref,
            ex_max=488,
            em_max=509,
            is_trusted=False,
            date_measured=date(2020, 1, 1),
        )
        m1.save(rebuild_cache=False)
        # Trusted but not from primary reference
        m2 = FM(
            fluorophore=state_with_ref,
            ex_max=500,
            em_max=520,
            is_trusted=True,
            date_measured=date(2010, 1, 1),
        )
        m2.save(rebuild_cache=False)

        state_with_ref.rebuild_attributes()
        state_with_ref.refresh_from_db()

        # Trusted measurement wins over primary reference
        assert state_with_ref.ex_max == 500
        assert state_with_ref.em_max == 520

    def test_date_priority_when_equal_trust_and_ref(self, state: State):
        """When trust and reference status are equal, more recent date wins."""
        m1 = FM(fluorophore=state, ex_max=488, em_max=509, date_measured=date(2020, 6, 1))
        m1.save(rebuild_cache=False)
        m2 = FM(fluorophore=state, ex_max=500, em_max=520, date_measured=date(2023, 6, 1))
        m2.save(rebuild_cache=False)

        state.rebuild_attributes()
        state.refresh_from_db()

        # More recent date wins
        assert state.ex_max == 500
        assert state.em_max == 520


class TestRebuildAttributesPinnedFields:
    """Test that pinned_source_map provides per-field override."""

    def test_pinned_field_overrides_all_priorities(self, state_with_ref: State):
        """A pinned field should override even trusted measurements."""
        primary_ref = state_with_ref.protein.primary_reference

        # Trusted measurement from primary reference - highest normal priority
        m_trusted = FM(
            fluorophore=state_with_ref,
            reference=primary_ref,
            ex_max=488,
            em_max=509,
            qy=0.67,
            is_trusted=True,
            date_measured=date(2020, 1, 1),
        )
        m_trusted.save(rebuild_cache=False)
        # Older, not trusted, not primary - lowest normal priority
        m_pinned = FM(
            fluorophore=state_with_ref,
            ex_max=500,
            em_max=520,
            qy=0.50,
            is_trusted=False,
            date_measured=date(2010, 1, 1),
        )
        m_pinned.save(rebuild_cache=False)

        # Pin qy to the low-priority measurement
        state_with_ref.pinned_source_map = {"qy": m_pinned.id}
        state_with_ref.save()

        state_with_ref.rebuild_attributes()
        state_with_ref.refresh_from_db()

        # ex_max and em_max should come from trusted measurement
        assert state_with_ref.ex_max == 488
        assert state_with_ref.em_max == 509
        # qy should come from the pinned measurement
        assert state_with_ref.qy == 0.50
        # source_map should reflect the pinned override
        assert state_with_ref.source_map["qy"] == m_pinned.id
        assert state_with_ref.source_map["ex_max"] == m_trusted.id

    def test_pinned_field_with_null_value_is_skipped(self, state: State):
        """If a pinned measurement has null for the pinned field, skip it."""
        m1 = FM(fluorophore=state, ex_max=488, qy=0.67, date_measured=date(2020, 1, 1))
        m1.save(rebuild_cache=False)
        m2 = FM(fluorophore=state, ex_max=500, qy=None, date_measured=date(2010, 1, 1))
        m2.save(rebuild_cache=False)

        # Pin qy to m2, but m2 has null qy
        state.pinned_source_map = {"qy": m2.id}
        state.save()

        state.rebuild_attributes()
        state.refresh_from_db()

        # qy should fall back to waterfall logic since pinned was null
        assert state.qy == 0.67
        assert state.source_map["qy"] == m1.id

    def test_pinned_nonexistent_measurement_is_ignored(self, state: State):
        """If pinned measurement doesn't exist, ignore and use waterfall."""
        m = FM(fluorophore=state, qy=0.67, date_measured=date(2020, 1, 1))
        m.save(rebuild_cache=False)

        # Pin to a non-existent measurement ID
        state.pinned_source_map = {"qy": 999999}
        state.save()

        state.rebuild_attributes()
        state.refresh_from_db()

        # Should fall back to waterfall logic
        assert state.qy == 0.67

    def test_pinned_measurement_from_other_fluorophore_is_ignored(self, protein: Protein):
        """A pinned measurement that belongs to a different fluorophore is ignored."""
        state1 = State.objects.create(protein=protein, name="state1")
        state2 = State.objects.create(protein=protein, name="state2")

        m1 = FM(fluorophore=state1, qy=0.67)
        m1.save(rebuild_cache=False)
        m_other = FM(fluorophore=state2, qy=0.50)
        m_other.save(rebuild_cache=False)

        # Try to pin state1's qy to state2's measurement
        state1.pinned_source_map = {"qy": m_other.id}
        state1.save()

        state1.rebuild_attributes()
        state1.refresh_from_db()

        # Should use state1's own measurement
        assert state1.qy == 0.67


class TestRebuildAttributesNoMeasurements:
    """Test behavior when there are no measurements."""

    def test_no_measurements_sets_null_values(self, state: State):
        """With no measurements, all values should be null."""
        state.rebuild_attributes()
        state.refresh_from_db()

        assert state.ex_max is None
        assert state.em_max is None
        assert state.qy is None
        assert state.source_map == {}


class TestRebuildAttributesAutoTrigger:
    """Test that save/delete on measurements auto-triggers rebuild."""

    def test_measurement_save_triggers_rebuild(self, state: State):
        """Creating a measurement should trigger rebuild_attributes."""
        # Create measurement - should auto-trigger rebuild (rebuild_cache=True by default)
        m = FM(fluorophore=state, ex_max=488, em_max=509)
        m.save()
        state.refresh_from_db()

        assert state.ex_max == 488
        assert state.em_max == 509

    def test_measurement_delete_triggers_rebuild(self, state: State):
        """Deleting a measurement should trigger rebuild_attributes."""
        m1 = FM(fluorophore=state, ex_max=488, em_max=509, date_measured=date(2023, 1, 1))
        m1.save()
        state.refresh_from_db()
        assert state.ex_max == 488

        m2 = FM(fluorophore=state, ex_max=500, em_max=520, date_measured=date(2020, 1, 1))
        m2.save()
        state.refresh_from_db()
        # Still 488 because more recent
        assert state.ex_max == 488

        # Delete the more recent measurement
        m1.delete()
        state.refresh_from_db()

        # Should now use the older measurement
        assert state.ex_max == 500
        assert state.em_max == 520

    def test_rebuild_cache_false_skips_rebuild(self, state: State):
        """rebuild_cache=False should skip the rebuild."""
        m = FM(fluorophore=state, ex_max=488)
        m.save(rebuild_cache=False)
        state.refresh_from_db()

        # Should not have been updated
        assert state.ex_max is None


class TestGetPrimaryReferenceId:
    """Test _get_primary_reference_id helper method."""

    def test_returns_none_for_no_primary_ref(self, state: State):
        """Should return None when protein has no primary reference."""
        assert state._get_primary_reference_id() is None

    def test_returns_primary_ref_id_for_protein(self, state_with_ref: State):
        """Should return primary_reference_id for protein state."""
        expected = state_with_ref.protein.primary_reference_id
        assert state_with_ref._get_primary_reference_id() == expected

    def test_returns_none_for_dye_without_ref(self, dyestate: DyeState):
        """Should return None when dye has no primary reference."""
        assert dyestate._get_primary_reference_id() is None

    def test_returns_primary_ref_id_for_dye(self, dyestate_with_ref: DyeState):
        """Should return primary_reference_id for dye state."""
        expected = dyestate_with_ref.dye.primary_reference_id
        assert dyestate_with_ref._get_primary_reference_id() == expected
