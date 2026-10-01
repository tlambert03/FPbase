"""The ajax views that write to the database take POSTs from logged-in users."""

from __future__ import annotations

from unittest.mock import patch

import pytest
from django.urls import reverse

from proteins.factories import ProteinFactory
from proteins.models import Protein
from references.factories import ReferenceFactory
from tests.test_users.factories import UserFactory

pytestmark = pytest.mark.django_db


@pytest.fixture
def protein() -> Protein:
    return ProteinFactory()


@pytest.fixture
def urls(protein: Protein) -> dict[str, str]:
    reference = ReferenceFactory()
    return {
        "filter_import": reverse("proteins:filter_import", args=("chroma",)),
        "add_reference": reverse("proteins:add_protein_reference", args=(protein.slug,)),
        "add_protein_excerpt": reverse("proteins:add_protein_excerpt", args=(protein.slug,)),
        "add_excerpt": reverse("references:add_excerpt", args=(reference.pk,)),
        "spectrum_preview": reverse("proteins:spectrum_preview"),
    }


VIEWS = [
    "filter_import",
    "add_reference",
    "add_protein_excerpt",
    "add_excerpt",
    "spectrum_preview",
]


@pytest.mark.parametrize("view", VIEWS)
def test_anonymous_post_is_sent_to_login(client, urls: dict[str, str], view: str):
    data = {"part": "ET525/50m", "reference_doi": "10.1038/nmeth.2413", "excerpt_content": "x"}
    with patch("proteins.views.spectra.add_filter_to_database") as add_filter:
        response = client.post(urls[view], data, HTTP_X_REQUESTED_WITH="XMLHttpRequest")
    assert response.status_code == 302
    assert response.url.startswith(reverse("account_login"))
    add_filter.assert_not_called()


@pytest.mark.parametrize("view", VIEWS)
def test_get_is_not_allowed(client, urls: dict[str, str], view: str):
    client.force_login(UserFactory())
    response = client.get(urls[view], HTTP_X_REQUESTED_WITH="XMLHttpRequest")
    assert response.status_code == 405


def test_approving_a_protein_takes_a_post(client, protein: Protein):
    Protein.objects.filter(id=protein.id).update(status=Protein.STATUS.pending)
    url = reverse("proteins:admin_approve_protein", args=(protein.slug,))
    client.force_login(UserFactory(is_staff=True))

    assert client.get(url).status_code == 405
    protein.refresh_from_db()
    assert protein.status == Protein.STATUS.pending

    assert client.post(url).status_code == 200
    protein.refresh_from_db()
    assert protein.status == Protein.STATUS.approved
