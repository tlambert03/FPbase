"""Object names are escaped in the links that admin pages build to related objects."""

from __future__ import annotations

import pytest
from django.urls import reverse
from django.utils.html import escape

from proteins.factories import MicroscopeFactory, OpticalConfigFactory, ProteinFactory
from proteins.models import ProteinCollection
from references.factories import ReferenceFactory
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db

NAME = "<b>name</b>"


def _assert_name_escaped(admin_client, url: str) -> None:
    content = admin_client.get(url).content.decode()
    assert escape(NAME) in content
    assert NAME not in content


def test_user_admin_links(admin_client):
    user = UserFactory()
    ProteinCollection.objects.create(name=NAME, owner=user)
    MicroscopeFactory(name=NAME, owner=user)
    _assert_name_escaped(admin_client, reverse("admin:users_user_change", args=(user.pk,)))


def test_reference_admin_links(admin_client):
    reference = ReferenceFactory()
    ProteinFactory(name=NAME, primary_reference=reference)
    url = reverse("admin:references_reference_change", args=(reference.pk,))
    _assert_name_escaped(admin_client, url)


def test_microscope_admin_links(admin_client):
    scope = OpticalConfigFactory(name=NAME).microscope
    _assert_name_escaped(
        admin_client, reverse("admin:proteins_microscope_change", args=(scope.pk,))
    )
