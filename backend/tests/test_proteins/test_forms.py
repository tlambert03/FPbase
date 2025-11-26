from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase

from proteins.forms import CollectionForm, ProteinForm, SpectrumForm, StateForm
from proteins.models import Protein, Spectrum, State

from ..test_users.factories import UserFactory

User = get_user_model()


class TestProteinForm(TestCase):
    def setUp(self):
        Protein.objects.get_or_create(
            name="Test Protein",
            seq="ARNDCEQGHILKMFPSTWYV",
            ipg_id="12345678",
            genbank="NC_000001.10",
            uniprot="P12345",
            pdb=["4HHB"],
        )

    def test_clean_proteinname_success(self):
        # Instantiate the form with a new protein
        form = ProteinForm({"name": "New Protein", "confirmation": True})
        # Run is_valid() to trigger the validation
        valid = form.is_valid()
        self.assertTrue(valid, "Form is not valid")

        # Run the actual clean_username method
        name = form.clean_name()
        self.assertEqual("New Protein", name)

    def test_clean_proteinname_exists(self):
        # Instantiate the form with existing protein name
        form = ProteinForm({"name": "Test Protein", "confirmation": True})
        # Run is_valid() to trigger the validation, which is going to fail
        # because the name is already taken
        valid = form.is_valid()
        self.assertFalse(valid)

        # The form.errors dict should contain a single error called 'name'
        self.assertTrue(len(form.errors) == 1)
        self.assertTrue("name" in form.errors)

    def test_clean_proteinseq_success(self):
        form = ProteinForm({"name": "New Protein", "seq": "ghilkmfpstwy varndceq", "confirmation": True})
        valid = form.is_valid()
        self.assertTrue(valid, "Form is not valid")
        seq = form.clean_seq()
        self.assertEqual("GHILKMFPSTWYVARNDCEQ", seq)

    def test_clean_proteinseq_exists(self):
        form = ProteinForm({"name": "New Protein", "seq": "ARNDCEQGHILKMFPSTWYV", "confirmation": True})
        valid = form.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form.errors) == 1)
        self.assertTrue("seq" in form.errors)

    def test_clean_proteinseq_invalid(self):
        form = ProteinForm({"name": "New Protein", "seq": "ARNDCEQGHILKMBZXFPSTWYV", "confirmation": True})
        valid = form.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form.errors) == 1)
        self.assertTrue("seq" in form.errors)

    def test_clean_refdoi_success(self):
        for doi in (
            "http://dx.doi.org/10.1038/nmeth.2413",
            "https://doi.org/10.1038/nmeth.2413",
            "10.1038/nmeth.2413",
        ):
            form = ProteinForm({"name": "New Protein", "reference_doi": doi, "confirmation": True})
            valid = form.is_valid()
            self.assertTrue(valid, "Form is not valid")
            self.assertEqual("10.1038/nmeth.2413", form.cleaned_data["reference_doi"])

    def test_clean_refdoi_failure(self):
        form = ProteinForm(
            {
                "name": "New Protein",
                "reference_doi": "30.1038/nmeth.2413",  # Invalid DOI
                "confirmation": True,
            }
        )
        valid = form.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form.errors) == 1)
        self.assertTrue("reference_doi" in form.errors)

    def test_ids_already_exist(self):
        form = ProteinForm(
            {
                "name": "New Protein",
                "ipg_id": "12345678",
                "genbank": "NC_000001.10",
                "uniprot": "P12345",
                "pdb": ["4HHB"],
                "confirmation": True,
            }
        )
        valid = form.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form.errors) == 4)
        self.assertTrue("ipg_id" in form.errors)
        self.assertTrue("genbank" in form.errors)
        self.assertTrue("uniprot" in form.errors)
        self.assertTrue("pdb" in form.errors)


class TestStateForm(TestCase):
    def setUp(self):
        self.t, _c = Protein.objects.get_or_create(name="Test Protein")
        State.objects.get_or_create(protein=self.t)

    def test_clean_state_success(self):
        form = StateForm({"name": "default", "ex_max": "488", "em_max": "525", "protein": self.t.id})
        valid = form.is_valid()
        self.assertTrue(valid, "Form is not valid")
        self.assertEqual("default", form.cleaned_data["name"])
        self.assertEqual(488, form.cleaned_data["ex_max"])
        self.assertEqual(525, form.cleaned_data["em_max"])

    def test_dark_state_success(self):
        form = StateForm({"name": "default", "is_dark": True, "protein": self.t.id})
        valid = form.is_valid()
        self.assertTrue(valid, "Form is not valid")

    def test_nondark_state_exemmax_required(self):
        form = StateForm({"name": "default", "is_dark": False, "protein": self.t.id})
        valid = form.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form.errors) == 2)
        self.assertTrue("ex_max" in form.errors)
        self.assertTrue("em_max" in form.errors)


class TestCollectionForm(TestCase):
    def setUp(self):
        self.p, _c = Protein.objects.get_or_create(name="Test Protein")
        self.userA = UserFactory(username="userA")
        self.userB = UserFactory(username="userB")

    def test_create_collection_success(self):
        form = CollectionForm(
            {
                "name": "New Collection",
                "description": "New collection description",
                "private": False,
            }
        )
        valid = form.is_valid()
        self.assertTrue(valid)

    def test_collection_clean_name_success(self):
        form = CollectionForm(
            {
                "name": "New Collection",
                "description": "New collection description",
                "private": False,
            },
            user=self.userA,
        )
        valid = form.is_valid()
        self.assertTrue(valid)

        # Run the actual clean_username method
        name = form.clean_name()
        self.assertEqual("New Collection", name)

    def test_collection_clean_name_failure(self):
        form = CollectionForm(
            {
                "name": "New Collection2",
                "description": "New collection description",
                "private": False,
            },
            user=self.userA,
        )
        valid = form.is_valid()
        self.assertTrue(valid)
        # attach user to instance
        obj = form.save(commit=False)
        obj.owner = self.userA
        obj.save()

        # userA cannot create another form with the same name
        form2 = CollectionForm(
            {
                "name": "New Collection2",
                "description": "New collection description",
                "private": False,
            },
            user=self.userA,
        )
        valid = form2.is_valid()
        self.assertFalse(valid)
        self.assertTrue(len(form2.errors) == 1)
        self.assertTrue("name" in form2.errors)

        # but let userB use the same Name
        form3 = CollectionForm(
            {
                "name": "New Collection2",
                "description": "New collection description",
                "private": False,
            },
            user=self.userB,
        )
        valid = form3.is_valid()
        self.assertTrue(valid)

        # let userA rename an existing form (update)
        form4 = CollectionForm(
            {
                "name": "new collection2",
                "description": "New collection description",
                "private": False,
            },
            user=self.userA,
            instance=obj,
        )
        valid = form4.is_valid()
        self.assertTrue(valid)


class TestSpectrumForm(TestCase):
    def setUp(self):
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

    def test_spectrum_form_manual_data_valid(self):
        """Test form validation with valid manual data and data_source=manual"""
        form_data = {
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
        form = SpectrumForm(data=form_data, files=None, user=self.user)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")

    def test_spectrum_form_manual_data_missing(self):
        """Test form validation with missing manual data when data_source=manual"""
        form_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",  # Empty manual data
            "data_source": "manual",
            "confirmation": True,
        }
        form = SpectrumForm(data=form_data, files=None, user=self.user)
        self.assertFalse(form.is_valid())
        self.assertIn("__all__", form.errors)
        self.assertIn("Please enter valid spectrum data", str(form.errors["__all__"]))

    def test_spectrum_form_file_data_valid(self):
        """Test form validation with valid file upload and data_source=file"""
        # Create a mock CSV file with consecutive wavelengths for step size = 1
        file_content = b"400,0.1\n401,0.2\n402,0.3\n403,0.5\n404,0.8\n405,1.0\n406,0.8\n407,0.5\n408,0.3\n409,0.1"
        uploaded_file = SimpleUploadedFile("spectrum.csv", file_content, content_type="text/csv")

        form_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",  # Empty manual data
            "data_source": "file",
            "confirmation": True,
        }
        files_data = {"file": uploaded_file}
        form = SpectrumForm(data=form_data, files=files_data, user=self.user)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")

    def test_spectrum_form_file_data_missing(self):
        """Test form validation with missing file when data_source=file"""
        form_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",
            "data_source": "file",
            "confirmation": True,
        }
        form = SpectrumForm(data=form_data, files={}, user=self.user)
        self.assertFalse(form.is_valid())
        self.assertIn("__all__", form.errors)
        self.assertIn("Please select a file to upload", str(form.errors["__all__"]))

    def test_spectrum_form_default_data_source(self):
        """Test that form defaults to file validation when no data_source is provided"""
        form_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": "",
            "confirmation": True,
            # No data_source provided - should default to "file"
        }
        form = SpectrumForm(data=form_data, files={}, user=self.user)
        self.assertFalse(form.is_valid())
        # Should show validation error since it defaults to file validation
        self.assertIn("__all__", form.errors)

    def test_spectrum_form_with_interpolation(self):
        """Test form validation with data requiring interpolation (step size > 1).

        This tests the fix for FPBASE-5G6 where NumPy float types from interpolation
        caused validation errors.
        """
        # Data with step size > 1 will trigger interpolation and NumPy floats
        form_data = {
            "category": Spectrum.PROTEIN,
            "subtype": Spectrum.EX,
            "owner_fluor": self.state.id,
            "data": ("[[400, 0.1], [405, 0.5], [410, 1.0], [415, 0.8], [420, 0.5], [425, 0.3], [430, 0.1]]"),
            "data_source": "manual",
            "confirmation": True,
        }
        form = SpectrumForm(data=form_data, files=None, user=self.user)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")
        # Ensure the form can be saved without NumPy type errors
        spectrum = form.save()
        self.assertIsNotNone(spectrum.id)
