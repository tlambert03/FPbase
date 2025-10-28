"""Pytest configuration for end-to-end tests.

This conftest.py is specifically for e2e tests that use Playwright and live_server.
These tests have different database cleanup requirements than unit tests.

Key differences from unit tests:
=============================
1. live_server runs Django in a separate thread, preventing transaction rollback
2. Database tables are TRUNCATED after each test (not rolled back)
3. Tests must be fully isolated - don't rely on session-scoped database fixtures

Why this matters:
================
- Unit tests use Django TestCase which rolls back transactions (fast, reliable)
- E2E tests with live_server CANNOT use transaction rollback (separate thread)
- pytest-django handles this by truncating all tables after each test
- Mixing these test types in one run can cause conflicts

References:
===========
- https://blog.tmk.name/2025/04/06/pytest-playwright-and-django/
- https://pytest-django.readthedocs.io/en/latest/database.html
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import warnings
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

import django.conf
import pytest
from allauth.account.models import EmailAddress
from django.contrib.auth import get_user_model
from django.contrib.sessions.backends.db import SessionStore
from webpack_loader import config, loaders, utils

if TYPE_CHECKING:
    from collections.abc import Iterator

    from django.contrib.auth.models import AbstractUser
    from playwright.sync_api import Browser, BrowserContext, ConsoleMessage, Page, ViewportSize

DEFAULT_TIMEOUT = int(os.getenv("DEFAULT_TIMEOUT", 8000))  # milliseconds
VIEWPORT_SIZE: ViewportSize = {"width": 1020, "height": 1200}

# Required for pytest-playwright fixtures to work with Django ORM in test setup.
# Playwright runs fixtures in an async event loop, triggering Django's async-unsafe
# protections. This is safe for tests because:
# 1) test fixtures run sequentially
# 2) live_server runs Django normally in a separate thread
# 3) no concurrent ORM access
# See: https://docs.djangoproject.com/en/5.2/topics/async/#envvar-DJANGO_ALLOW_ASYNC_UNSAFE
#      https://github.com/microsoft/playwright-pytest/issues/29
os.environ.setdefault("DJANGO_ALLOW_ASYNC_UNSAFE", "true")

# Mark all tests in this directory with transactional database access
# This is REQUIRED for live_server tests - it tells pytest-django to use
# table truncation instead of transaction rollback for cleanup
# See: https://pytest-django.readthedocs.io/en/latest/database.html#transactional-db
pytestmark = [
    pytest.mark.django_db(transaction=True),
]


def _frontend_assets_need_rebuild(stats_file) -> bool:
    """Check if frontend assets need to be rebuilt."""
    if not stats_file.is_file():
        return True

    assets = json.loads(stats_file.read_bytes())
    if assets.get("status") != "done" or not assets.get("chunks") or ("localhost" in assets.get("publicPath", "")):
        return True

    # Stats are valid - check if source files are newer
    stats_mtime = stats_file.stat().st_mtime
    frontend_src = Path(__file__).parent.parent.parent / "frontend" / "src"
    if any(
        f.stat().st_mtime > stats_mtime for f in frontend_src.rglob("*") if f.is_file() and not f.name.startswith(".")
    ):
        return True

    # Everything is up to date
    return False


@pytest.fixture(scope="module", autouse=True)
def _setup_frontend_assets() -> None:
    """Build webpack assets once per test module if needed.

    Checks if existing webpack stats represent a static build and if source files have changed.
    Rebuilds if missing, stale, dev server output detected, or source files are newer.

    Module-scoped to avoid rebuilding for every test (which would be slow).
    This is safe because webpack output is filesystem-based, not database.
    """
    stats_file = Path(django.conf.settings.WEBPACK_LOADER["DEFAULT"]["STATS_FILE"])

    # Need to build - either no stats, invalid stats, or source files changed
    if _frontend_assets_need_rebuild(stats_file):
        print("Building frontend assets for e2e tests...")
        subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)

    # UNDO the MockWebpackLoader used in normal unit tests in config.settings.test
    def _get_real_get_loader(config_name):
        return loaders.WebpackLoader(config_name, config.load_config(config_name))

    utils.get_loader = _get_real_get_loader


@pytest.fixture
def page(page: Page) -> Iterator[Page]:
    """Configure Playwright page fixture with FPbase defaults.

    Sets reasonable timeouts and viewport size for consistent test behavior.
    This wraps the default pytest-playwright page fixture.
    """
    # Set default timeout for actions (locator clicks, waits, etc)
    page.set_default_timeout(DEFAULT_TIMEOUT)
    # Set consistent viewport size for predictable rendering
    page.set_viewport_size(VIEWPORT_SIZE)
    with console_errors_raised(page):
        yield page


@contextmanager
def console_errors_raised(page: Page) -> Iterator[None]:
    """Context manager that collects console errors and warnings on a Playwright page."""
    messages: defaultdict[str, list[ConsoleMessage]] = defaultdict(list)
    ignore_messages = [
        "favicon.ico",
        "sentry",
        "WebGL",
        "[Report Only]",
    ]

    def on_console(msg: ConsoleMessage) -> None:
        for pattern in ignore_messages:
            if re.search(pattern, msg.text, re.IGNORECASE):
                return
        messages[msg.type].append(msg)

    page.on("console", on_console)

    yield

    if errors := messages.get("error"):
        error_messages = [f"  - {e.text} (at {e.location})" for e in errors]
        msg = "Console errors detected:\n" + "\n".join(error_messages)
        pytest.fail(msg)
    if _warnings := messages.get("warning"):
        warning_messages = [f"  - {w.text} (at {w.location})" for w in _warnings]
        msg = "Console warnings detected:\n" + "\n".join(warning_messages)
        warnings.warn(msg, stacklevel=2)


@pytest.fixture
def assert_no_console_errors(page: Page):
    """Fixture that collects console errors and asserts none occurred."""
    with console_errors_raised(page):
        yield


@pytest.fixture
def auth_user() -> AbstractUser:
    """Create authenticated user with verified email."""
    User = get_user_model()
    user = User.objects.create_user(username="testuser", email="test@example.com", password="testpass123")
    EmailAddress.objects.create(user=user, email=user.email, verified=True, primary=True)
    return user


@pytest.fixture
def auth_session_cookie(auth_user: AbstractUser) -> dict[str, str]:
    """Create Django session cookie for authenticated user."""
    session = SessionStore()
    session["_auth_user_id"] = str(auth_user.pk)
    session["_auth_user_backend"] = "django.contrib.auth.backends.ModelBackend"
    session["_auth_user_hash"] = auth_user.get_session_auth_hash()
    session.save()
    assert session.session_key is not None

    return {
        "name": "sessionid",
        "value": session.session_key,
        "domain": "localhost",
        "path": "/",
    }


@pytest.fixture
def auth_context(browser: Browser, auth_session_cookie: dict) -> Iterator[BrowserContext]:
    """Browser context with authenticated session."""
    context = browser.new_context(storage_state={"cookies": [auth_session_cookie]})
    yield context
    context.close()


@pytest.fixture
def auth_page(auth_context: BrowserContext) -> Iterator[Page]:
    """Page with authenticated session and console error checking."""
    page = auth_context.new_page()
    page.set_default_timeout(DEFAULT_TIMEOUT)
    page.set_viewport_size(VIEWPORT_SIZE)
    with console_errors_raised(page):
        yield page
    page.close()
