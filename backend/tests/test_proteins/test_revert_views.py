"""Tests for the staff-only revert endpoints."""

from __future__ import annotations

import pytest
import reversion
from django.test import Client
from django.urls import reverse
from reversion.models import Version

from proteins.models import Protein
from tests.test_users.factories import UserFactory


@pytest.fixture
def staff(db):
    return UserFactory(is_staff=True)


@pytest.fixture
def staff_client(staff):
    client = Client(enforce_csrf_checks=True)
    client.force_login(staff)
    return client


@pytest.fixture
def first_version(db) -> Version:
    with reversion.create_revision():
        protein = Protein.objects.create(name="Original", slug="original")
    with reversion.create_revision():
        protein.name = "Modified"
        protein.save()
    return Version.objects.get_for_object(protein).last()


def _urls(version: Version) -> list[str]:
    return [
        reverse("proteins:admin_revert_version", args=[version.id]),
        reverse("proteins:admin_revert_revision", args=[version.revision_id]),
    ]


@pytest.mark.parametrize("idx", [0, 1])
def test_revert_requires_post(staff_client, first_version, idx):
    assert staff_client.get(_urls(first_version)[idx]).status_code == 405
    assert Protein.objects.get(id=first_version.object_id).name == "Modified"


@pytest.mark.parametrize("idx", [0, 1])
def test_revert_requires_csrf_token(staff_client, first_version, idx):
    assert staff_client.post(_urls(first_version)[idx]).status_code == 403
    assert Protein.objects.get(id=first_version.object_id).name == "Modified"


@pytest.mark.parametrize("idx", [0, 1])
def test_revert_requires_staff(client, first_version, idx):
    client.force_login(UserFactory())
    client.post(_urls(first_version)[idx])
    assert Protein.objects.get(id=first_version.object_id).name == "Modified"


@pytest.mark.parametrize("idx", [0, 1])
def test_revert_records_new_revision(client, staff, first_version, idx):
    client.force_login(staff)
    protein = Protein.objects.get(id=first_version.object_id)
    n_versions = Version.objects.get_for_object(protein).count()

    assert client.post(_urls(first_version)[idx]).status_code == 200

    protein.refresh_from_db()
    assert protein.name == "Original"
    versions = Version.objects.get_for_object(protein)
    assert versions.count() == n_versions + 1
    assert versions.first().revision.user == staff
    assert "Reverted" in versions.first().revision.comment


def test_revert_unknown_version_404(client, staff):
    client.force_login(staff)
    url = reverse("proteins:admin_revert_version", args=[999999])
    assert client.post(url).status_code == 404
