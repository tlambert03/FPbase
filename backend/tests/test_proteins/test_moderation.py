"""Tests for what moderators see when users submit/approve things."""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest
from django.contrib.auth import get_user_model
from django.contrib.auth.models import Permission
from django.core import mail
from django.core.cache import cache
from django.urls import reverse

from proteins.factories import SpectrumFactory, StateFactory
from proteins.models import Spectrum
from proteins.models.spectrum import get_cached_spectra_info
from tests.test_proteins.test_views import INLINE_FORMSET

User = get_user_model()


@pytest.mark.django_db
def test_submission_email_describes_the_submission(client, settings) -> None:
    settings.MANAGERS = [("Manager", "manager@example.com")]
    user = User.objects.create_user(username="submitter", password="password")
    client.force_login(user)
    response = client.post(
        reverse("proteins:submit"),
        data={
            "name": "EmailedFP",
            "reference_doi": "10.1038/nmeth.2413",
            "states-0-name": "default",
            "states-0-ex_max": 488,
            "states-0-em_max": 525,
            "confirmation": True,
        }
        | INLINE_FORMSET,
    )
    assert response.status_code == 302

    (message,) = mail.outbox
    assert "User: submitter" in message.body
    assert "Protein: EmailedFP" in message.body
    assert response.url in message.body  # link to the protein page


@pytest.mark.django_db
def test_accepting_spectrum_updates_cached_spectra_list(client) -> None:
    moderator = User.objects.create_user(username="moderator", password="password")
    perms = Permission.objects.filter(codename__in=["change_spectrum", "delete_spectrum"])
    moderator.user_permissions.set(perms)
    client.force_login(moderator)

    spectrum = SpectrumFactory(
        owner_fluor=StateFactory(),
        category=Spectrum.PROTEIN,
        subtype=Spectrum.EX,
        status=Spectrum.STATUS.pending,
    )
    cache.clear()

    def cached_ids() -> set[int]:
        return {s["id"] for s in json.loads(get_cached_spectra_info())["data"]["spectra"]}

    assert spectrum.id not in cached_ids()  # pending: not listed (and the cache is now primed)

    with patch("proteins.views.spectra.uncache_protein_page") as uncache:
        response = client.post(
            reverse("proteins:pending_spectrum_action"),
            data={"spectrum_ids[]": [spectrum.id], "action": "accept"},
        )
    assert response.json()["success"]
    assert spectrum.id in cached_ids()
    # the owner protein's (view-cached) page must be refreshed too
    assert [c.args[0] for c in uncache.call_args_list] == [spectrum.owner_fluor.owner_slug]

    client.post(
        reverse("proteins:pending_spectrum_action"),
        data={"spectrum_ids[]": [spectrum.id], "action": "revert"},
    )
    assert spectrum.id not in cached_ids()
