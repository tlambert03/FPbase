"""Tests for object-level permissions on collection views."""

from __future__ import annotations

import pytest
from django.urls import reverse

from proteins.factories import ProteinFactory
from proteins.models import ProteinCollection
from tests.test_users.factories import UserFactory

AJAX = {"HTTP_X_REQUESTED_WITH": "XMLHttpRequest"}


@pytest.fixture
def owner(db):
    return UserFactory()


@pytest.fixture
def other_user(db):
    return UserFactory()


@pytest.fixture
def protein(db):
    return ProteinFactory()


@pytest.fixture
def private_collection(owner, protein):
    col = ProteinCollection.objects.create(name="secret", owner=owner, private=True)
    col.proteins.add(protein)
    return col


def _add(client, collection, protein):
    return client.post(
        reverse("proteins:add_to_collection"),
        {"collectionChoice": collection.id, "protein": protein.id},
        **AJAX,
    )


def test_add_to_collection_rejects_non_owner(client, other_user, private_collection):
    new_protein = ProteinFactory()
    client.force_login(other_user)
    response = _add(client, private_collection, new_protein)
    assert response.status_code == 403
    assert new_protein not in private_collection.proteins.all()


def test_add_to_collection_allows_owner_and_manager(client, owner, other_user):
    col = ProteinCollection.objects.create(name="shared", owner=owner, managers=[other_user.email])
    for user in (owner, other_user):
        new_protein = ProteinFactory()
        client.force_login(user)
        response = _add(client, col, new_protein)
        assert response.json() == {"status": "success"}
        assert new_protein in col.proteins.all()


@pytest.mark.parametrize("fmt", ["json", "csv", ""])
def test_private_collection_hidden_from_others(
    client, other_user, private_collection, protein, fmt
):
    url = reverse("proteins:collection-detail", args=[private_collection.id])
    for user in (None, other_user):
        if user:
            client.force_login(user)
        response = client.get(url, {"format": fmt} if fmt else {})
        content = b"".join(response) if response.streaming else response.content
        assert protein.name.encode() not in content
        assert b"This collection is not public" in content


@pytest.mark.parametrize("fmt", ["json", "csv", ""])
def test_private_collection_visible_to_owner(client, owner, private_collection, protein, fmt):
    client.force_login(owner)
    url = reverse("proteins:collection-detail", args=[private_collection.id])
    response = client.get(url, {"format": fmt} if fmt else {})
    content = b"".join(response) if response.streaming else response.content
    assert response.status_code == 200
    assert protein.name.encode() in content


def test_cannot_duplicate_private_collection(client, other_user, private_collection, protein):
    client.force_login(other_user)
    response = client.post(
        reverse("proteins:newcollection"),
        {"name": "copy", "dupcollection": private_collection.id},
    )
    assert response.status_code == 403
    assert not ProteinCollection.objects.filter(owner=other_user).exists()


def test_has_change_permission_does_not_mutate_managers(rf, owner):
    col = ProteinCollection.objects.create(name="mine", owner=owner)
    request = rf.get("/")
    request.user = owner
    assert col.has_change_permission(request)
    assert col.managers == []
