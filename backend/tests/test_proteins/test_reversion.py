"""Basic tests for django-reversion integration."""

import pytest
import reversion
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from reversion import is_registered
from reversion.models import Version

from proteins.models import Protein, State

User = get_user_model()


class TestReversionIntegration(TestCase):
    """Test that django-reversion 6.0 works with our models."""

    def test_models_registered_with_reversion(self):
        """Verify key models are registered with reversion."""
        assert is_registered(Protein), "Protein model should be registered"
        assert is_registered(State), "State model should be registered"

    def test_protein_update_creates_version(self):
        """Verify that updating a protein creates a version in the history."""
        with reversion.create_revision():
            protein = Protein.objects.create(name="TestProtein", slug="test-protein")

        # Get initial version count
        initial_versions = Version.objects.get_for_object(protein).count()

        # Update the protein in a revision context
        with reversion.create_revision():
            protein.name = "UpdatedProtein"
            protein.save()

        # Check that a new version was created
        updated_versions = Version.objects.get_for_object(protein).count()
        assert updated_versions > initial_versions, "Version should be created on update"

    def test_version_retrieval_and_revert(self):
        """Test retrieving version history and reverting to previous state."""
        with reversion.create_revision():
            protein = Protein.objects.create(name="OriginalName", slug="test-protein")
        original_name = protein.name

        # Update the protein in a revision context
        with reversion.create_revision():
            protein.name = "ModifiedName"
            protein.save()

        # Get version history
        versions = Version.objects.get_for_object(protein)
        assert versions.count() >= 1, "Should have at least one version"

        # Revert to first version
        first_version = versions.last()  # last() gets the oldest
        first_version.revision.revert()

        # Reload and check
        protein.refresh_from_db()
        assert protein.name == original_name, "Should revert to original name"


class TestReversionCompareAdmin(TestCase):
    """Test that django-reversion-compare works in admin."""

    def setUp(self):
        """Create a superuser for admin access."""
        self.user = User.objects.create_superuser(
            username="admin", email="admin@test.com", password="password"
        )
        self.client.login(username="admin", password="password")

    def test_admin_history_view_accessible(self):
        """Verify that admin history view loads without errors."""
        with reversion.create_revision():
            protein = Protein.objects.create(name="TestProtein", slug="test-protein")

        # Access the history view in admin
        url = reverse("admin:proteins_protein_history", args=[protein.pk])
        response = self.client.get(url)

        assert response.status_code == 200, "History view should be accessible"
        assert b"TestProtein" in response.content, "Protein name should appear in history"


@pytest.mark.django_db
def test_admin_approve_action_writes_snapshot(client) -> None:
    """Bulk-approving in the admin must leave an "approved" reversion snapshot."""
    admin = User.objects.create_superuser(username="admin", password="password")
    client.force_login(admin)
    with reversion.create_revision():
        protein = Protein.objects.create(name="PendingFP", status=Protein.STATUS.pending)
    before = protein.status_changed

    response = client.post(
        reverse("admin:proteins_protein_changelist"),
        {"action": "approve_protein", "_selected_action": [protein.pk]},
    )
    assert response.status_code == 302

    protein.refresh_from_db()
    assert protein.status == Protein.STATUS.approved
    assert protein.status_changed > before
    snapshot = Version.objects.get_for_object(protein).first()  # newest
    assert snapshot.field_dict["status"] == Protein.STATUS.approved
    assert snapshot.revision.user == admin
    assert "approved" in snapshot.revision.comment


@pytest.mark.django_db
def test_ajax_approve_writes_snapshot(client) -> None:
    staff = User.objects.create_user(username="staff", password="password", is_staff=True)
    client.force_login(staff)
    with reversion.create_revision():
        protein = Protein.objects.create(name="PendingFP", status=Protein.STATUS.pending)

    response = client.post(reverse("proteins:admin_approve_protein", args=[protein.slug]))
    assert response.status_code == 200

    protein.refresh_from_db()
    assert protein.status == Protein.STATUS.approved
    snapshot = Version.objects.get_for_object(protein).first()  # newest
    assert snapshot.field_dict["status"] == Protein.STATUS.approved
    assert snapshot.revision.user == staff
