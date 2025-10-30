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

import fcntl
import json
import os
import re
import subprocess
import warnings
from collections import defaultdict
from contextlib import contextmanager, suppress
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from allauth.account.models import EmailAddress
from django.contrib.auth import get_user_model
from django.contrib.sessions.backends.db import SessionStore
from playwright.sync_api import Page

# django-vite doesn't need loader imports

if TYPE_CHECKING:
    from collections.abc import Iterator

    from django.contrib.auth.models import AbstractUser
    from playwright.sync_api import Browser, BrowserContext, ConsoleMessage, Page, ViewportSize

# Register snapshot plugin to make assert_snapshot fixture available
pytest_plugins = ["tests_e2e.snapshot_plugin"]

DEFAULT_TIMEOUT = int(os.getenv("DEFAULT_TIMEOUT", 4000))  # milliseconds
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


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add custom command line options for e2e tests."""
    parser.addoption(
        "--visual-snapshots",
        action="store_true",
        default=False,
        help="Enable visual snapshot testing (disabled by default, never runs on CI)",
    )
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Update visual snapshots instead of comparing (for use with --visual-snapshots)",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Build frontend assets before test collection when running e2e tests.

    This hook runs very early - before test collection and before Django is imported.
    It ensures the manifest exists and is valid before any worker process starts.

    For pytest-xdist, this runs once in the main process before workers are spawned,
    ensuring all workers find a valid manifest when they start.
    """
    # Only run in main process (not in xdist workers)
    if hasattr(config, "workerinput"):
        return

    def _color_text(text: str, color_code: int) -> str:
        with suppress(AttributeError, OSError):
            # Check if stdout supports colors (has isatty and it returns True)
            if hasattr(sys.stdout, "isatty") and sys.stdout.isatty():
                return f"\033[{color_code}m{text}\033[0m"
        return text

    # Import here to avoid early Django import issues
    import django.conf

    manifest_file = Path(django.conf.settings.DJANGO_VITE["default"]["manifest_path"])
    lock_file = manifest_file.parent / ".build.lock"
    lock_file.parent.mkdir(parents=True, exist_ok=True)

    CYAN = 36
    RED = 31
    GREEN = 32

    # Use file locking to handle concurrent pytest runs
    with open(lock_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)

        if _frontend_assets_need_rebuild(manifest_file):
            print(_color_text("⚙️  Building frontend assets for e2e tests...", CYAN), flush=True)
            result = subprocess.run(
                ["pnpm", "--filter", "fpbase", "build"],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode != 0:
                print(_color_text(f"❌ Build failed with exit code {result.returncode}", RED), flush=True)
                print(f"STDOUT: {result.stdout}", flush=True)
                print(f"STDERR: {result.stderr}", flush=True)
                raise RuntimeError(f"Frontend build failed: {result.stderr}")
            print(_color_text("✅ Frontend build completed successfully", GREEN), flush=True)
        else:
            print(_color_text("✅ Frontend assets are up to date", GREEN), flush=True)


def _visual_snapshots_enabled(config: pytest.Config) -> bool:
    """Check if visual snapshots are enabled via CLI flag or environment variable."""
    return config.getoption("--visual-snapshots", False) or os.environ.get("VISUAL_SNAPSHOTS", "").lower() in (
        "1",
        "true",
        "yes",
    )


def _frontend_assets_need_rebuild(manifest_file) -> bool:
    """Check if frontend assets need to be rebuilt."""
    if not manifest_file.is_file():
        return True

    # Check if manifest is valid JSON
    try:
        manifest = json.loads(manifest_file.read_bytes())
        if not manifest:
            return True
    except (json.JSONDecodeError, ValueError):
        return True

    # Manifest is valid - check if source files are newer
    manifest_mtime = manifest_file.stat().st_mtime
    frontend_src = Path(__file__).parent.parent.parent / "frontend" / "src"
    if any(
        f.stat().st_mtime > manifest_mtime
        for f in frontend_src.rglob("*")
        if f.is_file() and not f.name.startswith(".")
    ):
        return True

    # Everything is up to date
    return False


# Frontend assets are built in pytest_configure hook (runs before Django import)
# No fixture needed here since assets are guaranteed to exist before workers start


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
