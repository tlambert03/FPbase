from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from proteins.models import Protein, State

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

        assert Protein.objects.count() == 0
        response = self.client.post(
            reverse("proteins:submit"),
            data={
                "name": "Test Protein",
                "reference_doi": "10.1038/nmeth.2413",
                "states-0-name": "default",
                "states-0-ex_max": 488,
                "states-0-em_max": 525,
                "confirmation": True,
            }
            | INLINE_FORMSET,
        )
        assert Protein.objects.count() == 1
        new_prot: Protein = Protein.objects.last()
        assert new_prot.name == "Test Protein"
        assert new_prot.primary_reference
        assert new_prot.primary_reference.doi == "10.1038/nmeth.2413"

        state: State = new_prot.default_state
        assert state.name == "default"
        assert state.ex_max == 488
        assert state.em_max == 525

        assert response.status_code == 302
        assert response.url == new_prot.get_absolute_url()
