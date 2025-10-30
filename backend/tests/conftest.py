import unittest.mock

import pytest


@pytest.fixture(scope="module", autouse=True)
def _mock_django_vite_for_unit_tests():
    """Mock django-vite asset loading for unit tests that don't need frontend assets.

    This fixture prevents django-vite from trying to load the manifest.json file
    when running unit tests without the frontend build. It returns empty strings
    for all asset requests, which is sufficient for tests that don't actually
    render templates or care about frontend assets.

    E2E tests (in tests_e2e/) skip this mock by having their own conftest.py
    that builds the frontend assets before tests run.
    """

    # Mock the generate_vite_asset method on the DjangoViteAssetLoader
    with unittest.mock.patch("django_vite.templatetags.django_vite.DjangoViteAssetLoader.instance") as mock_loader:
        mock_instance = unittest.mock.MagicMock()
        mock_instance.generate_vite_asset.return_value = ""
        mock_loader.return_value = mock_instance

        yield
