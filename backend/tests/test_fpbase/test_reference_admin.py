from __future__ import annotations

import pytest
from django.urls import reverse

from references.models import Reference

pytestmark = pytest.mark.django_db


def test_adding_a_reference_fetches_its_doi_info(admin_client) -> None:
    """The year etc. come from the DOI lookup, even with "refetch" left unchecked."""
    url = reverse("admin:references_reference_add")
    page = admin_client.get(url)
    data = {"doi": "10.1021/acs.biochem.3c00451"}
    for inline in page.context["inline_admin_formsets"]:
        prefix = inline.formset.prefix
        data |= {f"{prefix}-TOTAL_FORMS": "0", f"{prefix}-INITIAL_FORMS": "0"}

    response = admin_client.post(url, data)

    assert response.status_code == 302, response.context and response.context["errors"]
    reference = Reference.objects.get(doi=data["doi"])
    assert reference.year == 2024  # from the mocked DOI lookup
