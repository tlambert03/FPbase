from unittest.mock import patch

import pytest


@pytest.fixture(scope="session", autouse=True)
def _mock_django_vite_for_unit_tests():
    # Mock out the webpack loader to use a fake loader for unit tests
    # we do this here rather than config.settings.test to not interfere with e2e tests.
    from webpack_loader import loaders, utils

    with patch.object(utils, "get_loader", return_value=loaders.FakeWebpackLoader("DEFAULT", {})):
        yield
