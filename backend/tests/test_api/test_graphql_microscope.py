"""The microscope query fpbase-py sends (get_microscope): cost and output."""

from __future__ import annotations

import pytest
from django.db import connection
from django.test.utils import CaptureQueriesContext

from proteins.factories import (
    CameraFactory,
    LightFactory,
    MicroscopeFactory,
    OpticalConfigWithFiltersFactory,
)
from proteins.models import Spectrum

# as sent by fpbase-py's `get_microscope`
MICROSCOPE_QUERY = """
query getMicroscope($id: String!) {
    microscope(id: $id) {
        id
        opticalConfigs {
            name
            filters { path filter { id name spectrum { id subtype data } } }
            camera { id name spectrum { id subtype data } }
            light { id name spectrum { id subtype data } }
        }
    }
}
"""


def _query(client, scope_id: str) -> tuple[dict, int]:
    with CaptureQueriesContext(connection) as ctx:
        response = client.post(
            "/graphql/",
            {"query": MICROSCOPE_QUERY, "variables": {"id": scope_id}},
            content_type="application/json",
        )
    content = response.json()
    assert "errors" not in content, content["errors"]
    return content["data"]["microscope"], len(ctx.captured_queries)


def _scope_with_configs(n: int):
    scope = MicroscopeFactory()
    camera, light = CameraFactory(), LightFactory()
    for _ in range(n):
        OpticalConfigWithFiltersFactory(microscope=scope, camera=camera, light=light)
    return scope


@pytest.mark.django_db
def test_microscope_query_count_does_not_grow_with_configs(client):
    """Spectra come from prefetches, not a query per filter/camera/light."""
    _, small = _query(client, _scope_with_configs(1).id)
    _, large = _query(client, _scope_with_configs(5).id)
    assert large == small


@pytest.mark.django_db
def test_microscope_query_spectrum_data(client):
    """`data` is [[wavelength, value], ...], as with the old `[[Float]]` type."""
    scope = _scope_with_configs(1)
    data, _ = _query(client, scope.id)
    config = data["opticalConfigs"][0]
    owners = [fp["filter"] for fp in config["filters"]] + [config["camera"], config["light"]]
    for owner in owners:
        spectrum = Spectrum.objects.get(id=owner["spectrum"]["id"])
        assert owner["spectrum"]["data"] == [list(p) for p in spectrum.data]


@pytest.mark.django_db
def test_microscope_query_omits_unapproved_spectra(client):
    """A filter whose spectrum is pending review has no spectrum, as before."""
    scope = _scope_with_configs(1)
    placement = scope.optical_configs.first().filterplacement_set.first()
    Spectrum.objects.filter(owner_filter=placement.filter).update(status=Spectrum.STATUS.pending)
    data, _ = _query(client, scope.id)
    filters = {fp["filter"]["id"]: fp["filter"] for fp in data["opticalConfigs"][0]["filters"]}
    assert filters[str(placement.filter_id)]["spectrum"] is None
