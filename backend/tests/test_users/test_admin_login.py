"""The admin has no login form of its own: staff log in through the site's."""

from __future__ import annotations

import pytest
from django.urls import reverse

from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


def test_admin_login_redirects_to_site_login(client):
    response = client.get(reverse("admin:index"), follow=True)
    assert response.redirect_chain[-1][0].startswith(reverse("account_login"))

    credentials = {"username": "nobody", "password": "password"}
    response = client.post(reverse("admin:login"), credentials)
    assert response.status_code == 302
    assert response.url.startswith(reverse("account_login"))


def test_admin_login_refuses_non_staff(client):
    client.force_login(UserFactory())
    assert client.get(reverse("admin:login")).status_code == 403


def test_admin_is_available_to_logged_in_staff(client):
    client.force_login(UserFactory(is_staff=True))
    assert client.get(reverse("admin:index")).status_code == 200
