from __future__ import annotations

import pytest
from django.core import mail
from django.urls import reverse

from fpbase.users.models import User
from proteins.factories import MicroscopeFactory, ProteinFactory
from proteins.models import ProteinCollection
from references.factories import ReferenceFactory
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db
DELETE_URL = reverse("users:delete")


def _reauthenticate(client, user: User) -> None:
    client.force_login(user)
    resp = client.get(DELETE_URL)
    assert resp.status_code == 302
    assert resp.url.startswith(reverse("account_reauthenticate"))
    assert "Confirm Access" in client.get(resp.url).content.decode()
    resp = client.post(resp.url, {"password": "password"})
    assert resp.url == DELETE_URL


def test_change_username(client):
    user = UserFactory(username="old")
    client.force_login(user)
    resp = client.post(reverse("users:update"), {"username": "new", "name": "N"})
    assert resp.status_code == 302
    assert resp.url == reverse("users:detail", kwargs={"username": "new"})
    user.refresh_from_db()
    assert (user.username, user.name) == ("new", "N")
    assert len(mail.outbox) == 1
    assert mail.outbox[0].to == [user.email]
    assert '"old" to "new"' in mail.outbox[0].body


def test_change_name_only_sends_no_email(client):
    user = UserFactory(username="me")
    client.force_login(user)
    client.post(reverse("users:update"), {"username": "me", "name": "New Name"})
    assert not mail.outbox


def test_change_username_case_only(client):
    user = UserFactory(username="bob")
    client.force_login(user)
    client.post(reverse("users:update"), {"username": "Bob", "name": ""})
    user.refresh_from_db()
    assert user.username == "Bob"


def test_change_username_rejects_taken_name(client):
    UserFactory(username="taken")
    user = UserFactory(username="mine")
    client.force_login(user)
    resp = client.post(reverse("users:update"), {"username": "TAKEN", "name": ""})
    assert resp.status_code == 200
    assert resp.context["form"].errors["username"]
    user.refresh_from_db()
    assert user.username == "mine"


def test_delete_requires_login(client):
    resp = client.get(DELETE_URL)
    assert resp.status_code == 302
    assert resp.url.startswith(reverse("account_login"))


def test_delete_requires_reauthentication(client):
    user = UserFactory()
    client.force_login(user)
    resp = client.post(DELETE_URL, {"confirm": user.username})
    assert resp.url.startswith(reverse("account_reauthenticate"))
    assert User.objects.filter(pk=user.pk).exists()


def test_delete_requires_matching_username(client):
    user = UserFactory()
    _reauthenticate(client, user)
    for data in ({"confirm": "nope"}, {}):
        resp = client.post(DELETE_URL, data)
        assert resp.status_code == 200
        assert resp.context["form"].errors["confirm"]
    assert User.objects.filter(pk=user.pk).exists()


def test_delete_account(client):
    user = UserFactory()
    microscope = MicroscopeFactory(owner=user)
    collection = ProteinCollection.objects.create(name="mine", owner=user)
    protein = ProteinFactory(created_by=user)
    reference = ReferenceFactory(created_by=user, updated_by=user)

    _reauthenticate(client, user)
    profile = client.get(reverse("users:detail", kwargs={"username": user.username}))
    assert collection.name in profile.content.decode()
    assert "delete them first" in client.get(DELETE_URL).content.decode()

    email = user.email
    resp = client.post(DELETE_URL, {"confirm": user.username})
    assert resp.status_code == 302
    assert not User.objects.filter(pk=user.pk).exists()
    assert "_auth_user_id" not in client.session

    assert mail.outbox[-1].to == [email]
    assert "permanently deleted" in mail.outbox[-1].body

    # collections are deleted ...
    assert not ProteinCollection.objects.filter(pk=collection.pk).exists()
    # ... but everything else is kept, unattributed
    microscope.refresh_from_db()
    protein.refresh_from_db()
    reference.refresh_from_db()
    assert microscope.owner is None
    assert protein.created_by is None
    assert reference.created_by is None and reference.updated_by is None


def test_staff_cannot_delete_own_account(client):
    user = UserFactory(is_staff=True)
    _reauthenticate(client, user)
    resp = client.post(DELETE_URL, {"confirm": user.username})
    assert resp.url == reverse("users:detail", kwargs={"username": user.username})
    assert User.objects.filter(pk=user.pk).exists()
