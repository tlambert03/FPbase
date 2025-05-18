from test_plus.test import TestCase

from ..forms import CollectionForm, ProteinForm, StateForm
from ..models import Protein, State


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
        self.t, c = Protein.objects.get_or_create(name="Test Protein")
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
        self.p, c = Protein.objects.get_or_create(name="Test Protein")
        self.userA = self.make_user("userA", "userApassword")
        self.userB = self.make_user("userB", "userBpassword")

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
