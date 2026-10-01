"""Cached pages are stored per visitor."""

from __future__ import annotations

import pytest
import reversion
from django.core.cache import cache
from django.test import Client
from django.urls import reverse
from reversion.models import Version

from proteins.factories import ProteinFactory
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


@pytest.fixture(autouse=True)
def _clear_cache():
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def version_url() -> str:
    with reversion.create_revision():
        protein = ProteinFactory()
    version = Version.objects.get_for_object(protein).last()
    return f"{protein.get_absolute_url()}ver/{version.id}"


def _assert_cached_per_visitor(url: str) -> None:
    user = UserFactory()
    member, visitor = Client(), Client()
    member.force_login(user)
    for _ in range(2):  # (the first response sets the csrftoken cookie)
        assert user.email in member.get(url).content.decode()
    assert user.email not in visitor.get(url).content.decode()


def test_reference_list_cached_per_visitor():
    _assert_cached_per_visitor(reverse("reference-list"))


def test_protein_page_cached_per_visitor():
    _assert_cached_per_visitor(ProteinFactory().get_absolute_url())


def test_protein_version_page_cached_per_visitor(version_url: str):
    _assert_cached_per_visitor(version_url)
