"""End-to-end tests for the pending spectra moderation dashboard."""

from __future__ import annotations

from typing import TYPE_CHECKING

from django.contrib.auth.models import Permission
from django.urls import reverse
from playwright.sync_api import expect

from proteins.factories import SpectrumFactory, StateFactory
from proteins.models import Spectrum

if TYPE_CHECKING:
    from django.contrib.auth.models import AbstractUser
    from playwright.sync_api import Page
    from pytest_django.live_server_helper import LiveServer


def test_accept_and_undo_update_pending_count(
    live_server: LiveServer, auth_user: AbstractUser, auth_page: Page
) -> None:
    perms = Permission.objects.filter(codename__in=["change_spectrum", "delete_spectrum"])
    auth_user.user_permissions.set(perms)
    first, _second = (
        SpectrumFactory(owner_fluor=StateFactory(), category="p", subtype="ab", status="pending")
        for _ in range(2)
    )

    auth_page.goto(f"{live_server.url}{reverse('proteins:pending_spectra_dashboard')}")
    count = auth_page.locator("#pending-count")
    expect(count).to_have_text("2 pending")

    card = auth_page.locator(f'.spectrum-card[data-spectrum-id="{first.id}"]')
    card.get_by_role("button", name="Accept").click()

    toast = auth_page.locator(".notification-toast")
    expect(toast).to_have_count(1)
    expect(toast).to_contain_text("Accepted 1 spectrum(s)")
    expect(count).to_have_text("1 pending")
    expect(card).to_have_count(0)
    assert Spectrum.objects.all_objects().get(id=first.id).status == "approved"

    toast.get_by_role("button", name="Undo").click()
    expect(count).to_have_text("2 pending")
    expect(card).to_have_count(1)
    expect(auth_page.locator(".notification-toast.error")).to_have_count(0)
