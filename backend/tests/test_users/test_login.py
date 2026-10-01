"""Logging in goes through allauth's login form."""

from __future__ import annotations

import pytest
from django.urls import reverse

from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


def test_login_requires_verified_email(client):
    user = UserFactory()
    credentials = {"login": user.username, "password": "password"}
    response = client.post(reverse("account_login"), credentials)
    assert response.url == reverse("account_email_verification_sent")
    assert "_auth_user_id" not in client.session


def test_rest_framework_login_form_is_not_mounted(client):
    assert client.get("/api-auth/login/").status_code == 404
