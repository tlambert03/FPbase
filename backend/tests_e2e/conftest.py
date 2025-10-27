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
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import django.conf
import pytest

if TYPE_CHECKING:
    from playwright.sync_api import Page

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


@pytest.fixture(scope="module", autouse=True)
def _build_frontend_assets() -> None:
    """Build webpack assets once per test module if needed.

    Checks if existing webpack stats represent a static build.
    Rebuilds if missing, stale, or dev server output detected.

    Module-scoped to avoid rebuilding for every test (which would be slow).
    This is safe because webpack output is filesystem-based, not database.
    """
    stats_file = Path(django.conf.settings.WEBPACK_LOADER["DEFAULT"]["STATS_FILE"])
    if stats_file.is_file():
        assets = json.loads(stats_file.read_bytes())
        if (
            assets.get("status") == "done"
            and assets.get("chunks")
            and ("localhost" not in assets.get("publicPath", ""))
        ):
            return

    print("Building frontend assets for e2e tests...")
    subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)


@pytest.fixture
def page(page: Page) -> Page:
    """Configure Playwright page fixture with FPbase defaults.

    Sets reasonable timeouts and viewport size for consistent test behavior.
    This wraps the default pytest-playwright page fixture.
    """
    # Set default timeout for actions (locator clicks, waits, etc)
    page.set_default_timeout(10000)  # 10 seconds

    # Set consistent viewport size for predictable rendering
    page.set_viewport_size({"width": 1280, "height": 720})

    return page


@pytest.fixture
def assert_no_console_errors(page: Page):
    """Fixture that collects console errors and asserts none occurred.

    Usage:
        def test_something(page, assert_no_console_errors):
            page.goto(url)
            # Test interactions...
            # Errors are automatically checked at test end

    Filters out known acceptable errors:
    - Favicon 404s (browsers always request these)
    - Extension-related errors (from browser dev tools)

    Based on: https://playwright.dev/python/docs/api/class-page#page-event-console
    """
    # mapping of {level: [messages]}
    messages: defaultdict[str, list[str]] = defaultdict(list)
    ignore_messages = ["favicon.ico"]

    def on_console(msg):
        for pattern in ignore_messages:
            if re.search(pattern, msg.text, re.IGNORECASE):
                return
        messages[msg.type].append(msg)

    page.on("console", on_console)

    yield

    if errors := messages.get("error"):
        error_messages = [f"  - {e['text']} (at {e['location']})" for e in errors]
        msg = "Console errors detected:\n" + "\n".join(error_messages)
        pytest.fail(msg)
