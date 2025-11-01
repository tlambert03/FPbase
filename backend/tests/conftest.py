import pytest


@pytest.fixture(scope="module", autouse=True)
def _mock_django_vite_for_unit_tests():
    # Mock out the webpack loader to use a fake loader for unit tests
    # we do this here rather than config.settings.test to not interfere with e2e tests.
    from django.conf import settings
    from webpack_loader import utils

    settings.WEBPACK_LOADER["DEFAULT"]["LOADER_CLASS"] = "webpack_loader.loaders.FakeWebpackLoader"
    utils.get_loader.cache_clear()
