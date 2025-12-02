"""Test the external fpbase Python client against our GraphQL API.

This test ensures the fpbase PyPI package (https://pypi.org/project/fpbase/)
works correctly with our GraphQL API after schema changes.

The test handles a namespace conflict:
- `backend/fpbase/` is the Django app
- `fpbase` is also the external PyPI package

We use importlib to import the external package from its installed location.
"""

from __future__ import annotations

import importlib.util
import os
import site
import sys
from typing import TYPE_CHECKING

import pytest

from proteins.factories import DyeStateFactory, StateFactory
from proteins.models.fluorophore import FluorState

if TYPE_CHECKING:
    from pytest_django.live_server_helper import LiveServer


def _find_site_packages_fpbase():
    """Find the external fpbase package in site-packages.

    Since the Django app shadows the external package, we need to
    search sys.path manually for the site-packages version.
    """
    site_pkg = site.getsitepackages()[0]
    fpbase_init = f"{site_pkg}/fpbase/__init__.py"
    if os.path.exists(fpbase_init):
        return fpbase_init
    return None


def _import_fpbase_client():
    """Import the external fpbase package, working around namespace conflict.

    The Django app `backend/fpbase/` shadows the external `fpbase` package.
    We manually find and load the site-packages version.
    """
    # Find the external fpbase in site-packages
    fpbase_init = _find_site_packages_fpbase()
    if fpbase_init is None:
        pytest.skip("External fpbase package not installed. Install with: pip install fpbase")

    # Save current fpbase module if it exists (the Django app)
    django_fpbase = sys.modules.pop("fpbase", None)

    # Also need to clear submodules
    fpbase_submodules = {k: v for k, v in sys.modules.items() if k.startswith("fpbase.")}
    for k in fpbase_submodules:
        sys.modules.pop(k, None)

    try:
        # Create a spec from the site-packages path

        fpbase_dir = os.path.dirname(fpbase_init)
        spec = importlib.util.spec_from_file_location(
            "fpbase",
            fpbase_init,
            submodule_search_locations=[fpbase_dir],
        )

        if spec is None or spec.loader is None:
            pytest.skip("Could not create spec for external fpbase package")

        # Import the external package
        module = importlib.util.module_from_spec(spec)
        sys.modules["fpbase"] = module
        spec.loader.exec_module(module)

        return module
    finally:
        # Restore Django fpbase module and submodules if we saved them
        if django_fpbase is not None:
            sys.modules["fpbase"] = django_fpbase
        for k, v in fpbase_submodules.items():
            sys.modules[k] = v


@pytest.fixture
def fpbase_client(live_server: LiveServer):
    """Fixture that provides a configured fpbase client pointing to live_server."""
    _import_fpbase_client()

    # Import the internal _fetch module to configure the URL
    # Need to import it fresh since we may have restored the Django fpbase

    fpbase_init = _find_site_packages_fpbase()
    fpbase_dir = os.path.dirname(fpbase_init)
    fetch_path = os.path.join(fpbase_dir, "_fetch.py")

    spec = importlib.util.spec_from_file_location("fpbase._fetch", fetch_path)
    if spec is None or spec.loader is None:
        pytest.skip("Could not load fpbase._fetch module")

    fetch_module = importlib.util.module_from_spec(spec)
    sys.modules["fpbase._fetch"] = fetch_module
    spec.loader.exec_module(fetch_module)

    # Clear the singleton instance
    fetch_module.FPbaseClient._FPbaseClient__instance = None

    # Create a new client pointing at our live server
    client = fetch_module.FPbaseClient(base_url=f"{live_server.url}/graphql/")

    # Make it the singleton
    fetch_module.FPbaseClient._FPbaseClient__instance = client

    # Create a wrapper module-like object that has list_dyes, list_proteins, etc.
    class FPbaseClientWrapper:
        @staticmethod
        def list_dyes():
            return client.list_dyes()

        @staticmethod
        def list_proteins():
            return client.list_proteins()

    yield FPbaseClientWrapper()

    # Cleanup: clear the singleton
    fetch_module.FPbaseClient._FPbaseClient__instance = None
    sys.modules.pop("fpbase._fetch", None)


@pytest.mark.django_db(transaction=True)
def test_list_dyes_returns_dye_names(fpbase_client, live_server: LiveServer):
    """Test that fpbase.list_dyes() returns dye names from the API.

    This test validates that the GraphQL schema change from types.Dye
    to types.DyeState correctly returns dye data that the external
    fpbase client can parse.
    """

    # Create some dye states with distinct owner names
    # Use name="default" so the FluorState.label returns just the dye name
    dye_names = ["TestDyeA", "TestDyeB", "TestDyeC"]
    for name in dye_names:
        DyeStateFactory(dye__name=name, name=FluorState.DEFAULT_NAME)

    # Call list_dyes() through the external client
    result = fpbase_client.list_dyes()

    # Verify it's a list
    assert isinstance(result, list)

    # Verify our test dyes are in the results
    for name in dye_names:
        assert name in result, f"Expected {name!r} in list_dyes() result"


@pytest.mark.django_db(transaction=True)
def test_list_proteins_returns_protein_names(fpbase_client, live_server: LiveServer):
    """Test that fpbase.list_proteins() returns protein names from the API."""
    # Create some protein states
    protein_names = ["TestProtein1", "TestProtein2"]
    for name in protein_names:
        StateFactory(protein__name=name)

    # Call list_proteins() through the external client
    result = fpbase_client.list_proteins()

    # Verify it's a list
    assert isinstance(result, list)

    # Verify our test proteins are in the results
    for name in protein_names:
        assert name in result, f"Expected {name!r} in list_proteins() result"


@pytest.mark.django_db(transaction=True)
def test_list_dyes_no_none_values(fpbase_client, live_server: LiveServer):
    """Test that list_dyes() doesn't return None values (regression test).

    Before the fix, the GraphQL schema declared `dyes` as returning
    `List[types.Dye]` but the resolver returned `DyeState` objects.
    Graphene couldn't serialize them, returning None for each item.
    """
    # Create several dye states
    DyeStateFactory.create_batch(5)

    # Call list_dyes() - before the fix this would fail with:
    # TypeError: 'NoneType' object is not a mapping
    result = fpbase_client.list_dyes()

    # Should have results and no None values
    assert len(result) >= 5
    # The result is a list of strings (dye names), not dicts
    # If there were None values, the client code would have failed
    assert all(isinstance(name, str) for name in result)
