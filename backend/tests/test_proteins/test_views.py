import json
from typing import cast
from unittest.mock import patch

import pytest
from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from proteins.factories import (
    DyeStateFactory,
    MicroscopeFactory,
    OpticalConfigWithFiltersFactory,
    StateFactory,
)
from proteins.models import OcFluorEff, Protein, Spectrum, State
from proteins.tasks import calculate_scope_report

User = get_user_model()

INLINE_FORMSET = {
    "lineage-TOTAL_FORMS": 1,
    "lineage-INITIAL_FORMS": 0,
    "lineage-MIN_NUM_FORMS": 0,
    "lineage-MAX_NUM_FORMS": 1,
    "states-TOTAL_FORMS": 1,
    "states-INITIAL_FORMS": 0,
    "states-MIN_NUM_FORMS": 0,
    "states-MAX_NUM_FORMS": 1000,
}


class ProteinViewTests(TestCase):
    def setUp(self) -> None:
        self.admin_user = User.objects.create_superuser(
            username="admin", email="admin@example.com", password="password"
        )

    def test_protein_detail(self):
        """
        Test that the protein detail view returns a 200 response code
        """
        test_prot = Protein.objects.get_or_create(name="Test Protein")[0]
        response = self.client.get(test_prot.get_absolute_url())
        self.assertEqual(response.status_code, 200)

    def test_protein_submit(self):
        """
        Test that the protein detail view returns a 200 response code
        """
        self.client.login(username="admin", password="password")
        response = self.client.get(reverse("proteins:submit"))
        self.assertEqual(response.status_code, 200)

        name = "Protein ERMCOFSD"
        initial_count = Protein.objects.count()
        response = self.client.post(
            reverse("proteins:submit"),
            data={
                "name": name,
                "reference_doi": "10.1038/nmeth.2413",
                "states-0-name": "default",
                "states-0-ex_max": 488,
                "states-0-em_max": 525,
                "confirmation": True,
            }
            | INLINE_FORMSET,
        )
        assert response.status_code == 302

        assert Protein.objects.count() == initial_count + 1
        new_prot = cast("Protein", Protein.objects.get(name=name))
        assert response.url == new_prot.get_absolute_url()

        assert new_prot.name == name
        assert new_prot.primary_reference
        assert new_prot.primary_reference.doi == "10.1038/nmeth.2413"

        state = new_prot.default_state
        assert isinstance(state, State)
        assert state.name == "default"
        assert state.ex_max == 488
        assert state.em_max == 525


class SpectrumPreviewViewTests(TestCase):
    def setUp(self) -> None:
        self.user = User.objects.create_user(username="testuser", password="testpass")
        self.protein = Protein.objects.create(
            name="Test Protein",
            seq="ARNDCEQGHILKMFPSTWYV",
            ipg_id="12345678",
        )
        self.state = State.objects.create(
            name="default",
            protein=self.protein,
            ex_max=488,
            em_max=525,
        )
        self.preview_url = reverse("proteins:spectrum_preview")

    def test_spectrum_preview_requires_post(self):
        """Test that spectrum preview endpoint requires POST method"""
        self.client.login(username="testuser", password="testpass")
        response = self.client.get(self.preview_url)
        self.assertEqual(response.status_code, 405)
        data = json.loads(response.content)
        self.assertEqual(data["error"], "POST required")

    def test_spectrum_preview_manual_data_success(self):
        """Test successful spectrum preview with manual data"""
        self.client.login(username="testuser", password="testpass")

        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": (
                "[[400, 0.1], [401, 0.2], [402, 0.3], [403, 0.5], [404, 0.8], "
                "[405, 1.0], [406, 0.8], [407, 0.5], [408, 0.3], [409, 0.1]]"
            ),
            "data_source": "manual",
            "confirmation": True,
        }

        response = self.client.post(self.preview_url, data=post_data)
        self.assertEqual(response.status_code, 200)

        data = json.loads(response.content)
        self.assertTrue(data["success"])
        self.assertIn("preview", data)
        self.assertIn("svg", data["preview"])
        self.assertIn("peak_wave", data["preview"])
        self.assertIn("data_points", data["preview"])
        self.assertFalse(data["preview"]["file_was_processed"])

    def test_spectrum_preview_file_upload_success(self):
        """Test successful spectrum preview with file upload"""
        self.client.login(username="testuser", password="testpass")

        # Create a mock CSV file with consecutive wavelengths
        file_content = (
            b"400,0.1\n401,0.2\n402,0.3\n403,0.5\n404,0.8\n"
            b"405,1.0\n406,0.8\n407,0.5\n408,0.3\n409,0.1"
        )
        uploaded_file = SimpleUploadedFile("spectrum.csv", file_content, content_type="text/csv")

        # Use multipart form data for file upload
        response = self.client.post(
            self.preview_url,
            data={
                "category": Spectrum.PROTEIN,
                "subtype": Spectrum.EX,
                "owner_fluor": self.state.id,
                "data": "",
                "data_source": "file",
                "confirmation": True,
                "file": uploaded_file,
            },
            format="multipart",
        )
        self.assertEqual(response.status_code, 200)

        data = json.loads(response.content)
        self.assertTrue(data["success"])
        self.assertIn("preview", data)
        self.assertIn("svg", data["preview"])
        self.assertTrue(data["preview"]["file_was_processed"])

    def test_spectrum_preview_validation_failure_manual(self):
        """Test spectrum preview with invalid manual data"""
        self.client.login(username="testuser", password="testpass")

        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",  # Empty manual data
            "data_source": "manual",
            "confirmation": True,
        }

        response = self.client.post(self.preview_url, data=post_data)
        self.assertEqual(response.status_code, 400)

        data = json.loads(response.content)
        self.assertEqual(data["error"], "Form validation failed. Please check your input data.")
        self.assertIn("form_errors", data)
        # Error should be in non-field errors (__all__)
        self.assertTrue("__all__" in data["form_errors"] or "data" in data["form_errors"])

    def test_spectrum_preview_validation_failure_file(self):
        """Test spectrum preview with missing file upload"""
        self.client.login(username="testuser", password="testpass")

        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",
            "data_source": "file",
            "confirmation": True,
        }

        response = self.client.post(self.preview_url, data=post_data)
        self.assertEqual(response.status_code, 400)

        data = json.loads(response.content)
        self.assertEqual(data["error"], "Form validation failed. Please check your input data.")
        self.assertIn("form_errors", data)
        # Error should be in non-field errors (__all__)
        self.assertTrue("__all__" in data["form_errors"] or "file" in data["form_errors"])

    def test_spectrum_preview_invalid_spectrum_data(self):
        """Test spectrum preview with invalid spectrum data format"""
        self.client.login(username="testuser", password="testpass")

        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "invalid data format",  # Invalid spectrum data
            "data_source": "manual",
            "confirmation": True,
        }

        response = self.client.post(self.preview_url, data=post_data)
        self.assertEqual(response.status_code, 400)

        data = json.loads(response.content)
        # Should be either form validation error or data processing error
        self.assertIn("error", data)

    def test_spectrum_preview_requires_authentication(self):
        """Test that spectrum preview fails for anonymous users"""
        # Don't log in - test anonymous user
        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": (
                "[[400, 0.1], [401, 0.2], [402, 0.3], [403, 0.5], [404, 0.8], "
                "[405, 1.0], [406, 0.8], [407, 0.5], [408, 0.3], [409, 0.1]]"
            ),
            "data_source": "manual",
            "confirmation": True,
        }

        response = self.client.post(self.preview_url, data=post_data)
        # Should fail with 500 error because anonymous user can't be assigned to created_by
        self.assertEqual(response.status_code, 500)

        data = json.loads(response.content)
        self.assertIn("error", data)

    def test_spectrum_preview_data_source_defaults_to_file(self):
        """Test that data_source defaults to 'file' when not provided"""
        self.client.login(username="testuser", password="testpass")

        post_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",
            "confirmation": True,
            # No data_source provided
        }

        response = self.client.post(self.preview_url, data=post_data)
        self.assertEqual(response.status_code, 400)

        data = json.loads(response.content)
        self.assertIn("form_errors", data)
        # Should show validation error since it defaults to file validation
        self.assertTrue("__all__" in data["form_errors"] or "file" in data["form_errors"])


@pytest.mark.django_db
class TestScopeReportJson:
    """Test scope_report_json view that serves microscope efficiency report data.

    These tests exercise the query paths through the Fluorophore MTI structure,
    ensuring both State (protein) and DyeState entries are handled correctly.
    """

    def test_scope_report_json_returns_valid_response(self, client):
        """Test that scope_report_json returns valid JSON structure."""
        microscope = MicroscopeFactory()
        OpticalConfigWithFiltersFactory(microscope=microscope)

        url = reverse("proteins:scope_report_json", args=[microscope.id])
        response = client.get(url)

        assert response.status_code == 200
        data = response.json()

        # Check expected top-level structure
        assert "microscope" in data
        assert "report" in data
        assert "fluors" in data
        assert data["microscope"] == microscope.id

    def test_scope_report_json_includes_protein_states(self, client):
        """Test that scope_report_json correctly includes protein States."""
        microscope = MicroscopeFactory()
        oc = OpticalConfigWithFiltersFactory(microscope=microscope)
        state = StateFactory()

        # Create OcFluorEff record using bulk_create to bypass save() which
        # would call update_effs() and overwrite our test values
        OcFluorEff.objects.bulk_create(
            [
                OcFluorEff(
                    oc=oc,
                    fluor=state,
                    fluor_name=str(state),
                    ex_eff=0.8,
                    em_eff=0.7,
                    brightness=50.0,
                )
            ]
        )

        url = reverse("proteins:scope_report_json", args=[microscope.id])
        response = client.get(url)

        assert response.status_code == 200
        data = response.json()

        # State should appear in report
        assert len(data["report"]) > 0
        config_data = data["report"][0]
        assert "values" in config_data

        # Find the state in the values
        state_entries = [v for v in config_data["values"] if v["fluor_slug"] == state.slug]
        assert len(state_entries) == 1, f"State {state.slug} should appear in report"

        entry = state_entries[0]
        assert entry["shape"] == "circle", "Protein states should have circle shape"
        assert entry["ex_eff"] == 0.8
        assert entry["em_eff"] == 0.7

        # State should appear in fluors dict
        assert state.slug in data["fluors"]
        fluor_data = data["fluors"][state.slug]
        assert fluor_data["type"] == "p", "Protein fluorophores should have type 'p'"

    def test_scope_report_json_includes_dye_states(self, client):
        """Test that scope_report_json correctly includes DyeStates."""
        microscope = MicroscopeFactory()
        oc = OpticalConfigWithFiltersFactory(microscope=microscope)
        dye_state = DyeStateFactory()

        # Create OcFluorEff record using bulk_create to bypass save()
        OcFluorEff.objects.bulk_create(
            [
                OcFluorEff(
                    oc=oc,
                    fluor=dye_state,
                    fluor_name=str(dye_state),
                    ex_eff=0.9,
                    em_eff=0.6,
                    brightness=40.0,
                )
            ]
        )

        url = reverse("proteins:scope_report_json", args=[microscope.id])
        response = client.get(url)

        assert response.status_code == 200
        data = response.json()

        # DyeState should appear in report
        assert len(data["report"]) > 0
        config_data = data["report"][0]

        # Find the dye state in the values
        dye_entries = [v for v in config_data["values"] if v["fluor_slug"] == dye_state.slug]
        assert len(dye_entries) == 1, f"DyeState {dye_state.slug} should appear in report"

        entry = dye_entries[0]
        assert entry["shape"] == "square", "Dye states should have square shape"
        assert entry["ex_eff"] == 0.9
        assert entry["em_eff"] == 0.6

        # DyeState should appear in fluors dict
        assert dye_state.slug in data["fluors"]
        fluor_data = data["fluors"][dye_state.slug]
        assert fluor_data["type"] == "d", "Dye fluorophores should have type 'd'"

    def test_scope_report_json_with_both_states_and_dyes(self, client):
        """Test that scope_report_json handles mixed State and DyeState data.

        This is the key test for validating the migration's handling of
        the new Fluorophore MTI structure with both entity types.
        """
        microscope = MicroscopeFactory()
        oc = OpticalConfigWithFiltersFactory(microscope=microscope)
        state = StateFactory()
        dye_state = DyeStateFactory()

        # Create OcFluorEff records for both using bulk_create
        OcFluorEff.objects.bulk_create(
            [
                OcFluorEff(oc=oc, fluor=state, fluor_name=str(state), ex_eff=0.8, em_eff=0.7),
                OcFluorEff(
                    oc=oc, fluor=dye_state, fluor_name=str(dye_state), ex_eff=0.9, em_eff=0.6
                ),
            ]
        )

        url = reverse("proteins:scope_report_json", args=[microscope.id])
        response = client.get(url)

        assert response.status_code == 200
        data = response.json()

        # Both should appear in report
        config_data = data["report"][0]
        slugs_in_report = [v["fluor_slug"] for v in config_data["values"]]

        assert state.slug in slugs_in_report, "State should appear in report"
        assert dye_state.slug in slugs_in_report, "DyeState should appear in report"

        # Both should be in fluors dict with correct types
        assert data["fluors"][state.slug]["type"] == "p"
        assert data["fluors"][dye_state.slug]["type"] == "d"

    def test_scope_report_json_full_workflow_integration(self, client):
        """Integration test: run calculate_scope_report then verify JSON output.

        This tests the complete workflow from task execution to JSON response,
        validating that the migration preserves the end-to-end behavior.
        """
        microscope = MicroscopeFactory()
        OpticalConfigWithFiltersFactory(microscope=microscope)
        state = StateFactory()
        dye_state = DyeStateFactory()

        # Verify both have spectra
        assert state.has_spectra()
        assert dye_state.has_spectra()

        # Run the task to create OcFluorEff records
        with patch.object(calculate_scope_report, "update_state"):
            calculate_scope_report.run(scope_id=microscope.id)

        # Verify OcFluorEff records were created
        assert OcFluorEff.objects.count() >= 2

        # Now test the JSON endpoint
        url = reverse("proteins:scope_report_json", args=[microscope.id])
        response = client.get(url)

        assert response.status_code == 200
        data = response.json()

        # Verify structure
        assert "report" in data
        assert "fluors" in data

        # Verify fluors dict contains our fluorophores
        assert state.slug in data["fluors"], f"State {state.slug} missing from fluors"
        assert dye_state.slug in data["fluors"], f"DyeState {dye_state.slug} missing from fluors"


class TestProteinSubmitMultiState(TestCase):
    """Regression tests for FPBASE-6H7: protein submission with multiple states."""

    def setUp(self) -> None:
        self.admin_user = User.objects.create_superuser(
            username="admin", email="admin@example.com", password="password"
        )

    def test_protein_submit_with_two_states(self):
        """Test protein submission with multiple states doesn't raise ValueError.

        Regression test for FPBASE-6H7: When a protein has 2+ states but no transitions,
        the check_switch_type function is called and previously crashed with ValueError
        due to incorrect dict(Protein.SwitchingChoices) usage.
        """
        self.client.login(username="admin", password="password")
        initial_count = Protein.objects.count()

        # Submit protein with 2 states (no transitions) - triggers check_switch_type
        # with suggested switch_type='o' (OTHER)
        response = self.client.post(
            reverse("proteins:submit"),
            data={
                "name": "MultiStateProtein",
                "reference_doi": "10.1038/nmeth.2413",
                # Two states, no transitions = triggers 'OTHER' suggestion
                "states-0-name": "on",
                "states-0-ex_max": 488,
                "states-0-em_max": 525,
                "states-1-name": "off",
                "states-1-ex_max": 400,
                "states-1-em_max": 450,
                "confirmation": True,
                "lineage-TOTAL_FORMS": 1,
                "lineage-INITIAL_FORMS": 0,
                "lineage-MIN_NUM_FORMS": 0,
                "lineage-MAX_NUM_FORMS": 1,
                "states-TOTAL_FORMS": 2,
                "states-INITIAL_FORMS": 0,
                "states-MIN_NUM_FORMS": 0,
                "states-MAX_NUM_FORMS": 1000,
            },
        )

        # Should redirect to protein detail (success), not error
        assert response.status_code == 302
        assert Protein.objects.count() == initial_count + 1

        new_prot = cast("Protein", Protein.objects.get(name="MultiStateProtein"))
        assert response.url == new_prot.get_absolute_url()
        assert new_prot.states.count() == 2
