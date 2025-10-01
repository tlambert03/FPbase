import json
from typing import cast

from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from proteins.models import Protein, Spectrum, State

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

        state: State = new_prot.default_state
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
            "owner_state": self.state.id,
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
        file_content = b"400,0.1\n401,0.2\n402,0.3\n403,0.5\n404,0.8\n405,1.0\n406,0.8\n407,0.5\n408,0.3\n409,0.1"
        uploaded_file = SimpleUploadedFile("spectrum.csv", file_content, content_type="text/csv")

        # Use multipart form data for file upload
        response = self.client.post(
            self.preview_url,
            data={
                "category": Spectrum.PROTEIN,
                "subtype": Spectrum.EX,
                "owner_state": self.state.id,
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
            "owner_state": self.state.id,
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
            "owner_state": self.state.id,
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
            "owner_state": self.state.id,
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
            "owner_state": self.state.id,
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
            "owner_state": self.state.id,
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
