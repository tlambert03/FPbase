from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from django.test import Client
    from pytest_django.fixtures import SettingsWrapper


@pytest.mark.django_db
def test_ga_snippet_requires_setting(client: Client, settings: SettingsWrapper) -> None:
    """Google Analytics must only load when GOOGLE_ANALYTICS_ID is set (i.e. production)."""
    assert "googletagmanager" not in client.get("/").content.decode()

    settings.GOOGLE_ANALYTICS_ID = "G-TEST"
    assert "gtag/js?id=G-TEST" in client.get("/").content.decode()
