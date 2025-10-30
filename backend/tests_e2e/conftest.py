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
import shutil
import subprocess
import sys
import warnings
from collections import defaultdict
from contextlib import contextmanager, suppress
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import django.conf
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

    CYAN = 36
    RED = 31
    GREEN = 32

    def _color_text(text: str, color_code: int) -> str:
        with suppress(AttributeError, OSError):
            # Check if stdout supports colors (has isatty and it returns True)
            if hasattr(sys.stdout, "isatty") and sys.stdout.isatty():
                return f"\033[{color_code}m{text}\033[0m"
        return text

    manifest_file = Path(django.conf.settings.DJANGO_VITE["default"]["manifest_path"])
    lock_file = manifest_file.parent / ".build.lock"
    lock_file.parent.mkdir(parents=True, exist_ok=True)

    # Use file locking to handle concurrent pytest runs
    with open(lock_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)

        if _frontend_assets_need_rebuild(manifest_file):
            print(_color_text("🔨 Building frontend assets for e2e tests...", CYAN), flush=True)
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

            # Clear django-vite's cached manifest so it reloads the new one
            # This is necessary because Django may have already loaded the old manifest
            import django_vite.core.asset_loader

            django_vite.core.asset_loader.DjangoViteAssetLoader._instance = None
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


# Add a data store for computed paths
class SnapshotPaths:
    snapshots_path: Path | None = None
    failures_path: Path | None = None


@pytest.fixture(scope="session", autouse=True)
def _cleanup_snapshot_failures(pytestconfig: pytest.Config):
    """
    Clean up snapshot failures directory once at the beginning of test session.

    The snapshot storage path is relative to each test folder, modeling after the React snapshot locations
    """

    root_dir = Path(pytestconfig.rootdir)

    # Compute paths once
    SnapshotPaths.snapshots_path = root_dir / "__snapshots__"

    SnapshotPaths.failures_path = root_dir / "snapshot_failures"

    # Clean up the entire failures directory at session start so past failures don't clutter the result
    if SnapshotPaths.failures_path.exists():
        shutil.rmtree(SnapshotPaths.failures_path, ignore_errors=True)

    # Create the directory to ensure it exists
    with suppress(FileExistsError):
        SnapshotPaths.failures_path.mkdir(parents=True, exist_ok=True)

    yield


@pytest.fixture
def assert_snapshot(pytestconfig: pytest.Config, request: pytest.FixtureRequest) -> Callable:
    if not _visual_snapshots_enabled(pytestconfig):

        def noop(*args: Any, **kwargs: Any) -> None:
            pass

        noop.NOOP = True  # type: ignore[attr-defined]s
        return noop

    from io import BytesIO

    import pytest
    from PIL import Image
    from pixelmatch.contrib.PIL import pixelmatch

    test_function_name = request.node.name
    SNAPSHOT_MESSAGE_PREFIX = "[playwright-visual-snapshot]"

    test_name_without_params = test_function_name.split("[", 1)[0]
    test_name = f"{test_function_name}[{sys.platform!s}]"

    current_test_file_path = Path(request.node.fspath)
    current_test_file_path.parent.resolve()

    # Use global paths if available, otherwise calculate per test
    snapshots_path = SnapshotPaths.snapshots_path
    assert snapshots_path

    snapshot_failures_path = SnapshotPaths.failures_path
    assert snapshot_failures_path

    # we know this exists because of the default value on ini
    global_snapshot_threshold = 0.1

    mask_selectors = []
    update_snapshot = pytestconfig.getoption("--update-snapshots", False)

    # for automatically naming multiple assertions
    counter = 0
    # Collection to store failures
    failures = []

    def _create_locators_from_selectors(page: Page, selectors: list[str]):
        """
        Convert a list of CSS selector strings to locator objects
        """
        return [page.locator(selector) for selector in selectors]

    def compare(
        img_or_page: bytes | Any,
        *,
        threshold: float | None = None,
        name=None,
        fail_fast=False,
        mask_elements: list[str] | None = None,
    ) -> None:
        nonlocal counter

        if not name:
            if counter > 0:
                name = f"{test_name}_{counter}.png"
            else:
                name = f"{test_name}.png"

        # Use global threshold if no local threshold provided
        if not threshold:
            threshold = global_snapshot_threshold

        # If page reference is passed, use screenshot
        if isinstance(img_or_page, Page):
            # Combine configured mask elements with any provided in the function call
            all_mask_selectors = list(mask_selectors)
            if mask_elements:
                all_mask_selectors.extend(mask_elements)

            # Convert selectors to locators
            masks = _create_locators_from_selectors(img_or_page, all_mask_selectors) if all_mask_selectors else []

            img = img_or_page.screenshot(
                animations="disabled",
                type="png",
                mask=masks,
                # TODO only for jpeg
                # quality=100,
            )
        else:
            img = img_or_page

        # test file without the extension
        test_file_name_without_extension = current_test_file_path.stem

        # Created a nested folder to store screenshots: snapshot/test_file_name/test_name/
        test_file_snapshot_dir = snapshots_path / test_file_name_without_extension / test_name_without_params
        test_file_snapshot_dir.mkdir(parents=True, exist_ok=True)

        screenshot_file = test_file_snapshot_dir / name

        # Create a dir where all snapshot test failures will go
        # ex: snapshot_failures/test_file_name/test_name
        failure_results_dir = snapshot_failures_path / test_file_name_without_extension / test_name

        # increment counter before any failures are recorded
        counter += 1

        if update_snapshot:
            screenshot_file.write_bytes(img)
            failures.append(f"{SNAPSHOT_MESSAGE_PREFIX} Snapshots updated. Please review images. {screenshot_file}")
            return

        if not screenshot_file.exists():
            screenshot_file.write_bytes(img)
            failures.append(
                f"{SNAPSHOT_MESSAGE_PREFIX} New snapshot(s) created. Please review images. {screenshot_file}"
            )
            return

        img_a = Image.open(BytesIO(img))
        img_b = Image.open(screenshot_file)
        img_diff = Image.new("RGBA", img_a.size)
        mismatch = pixelmatch(img_a, img_b, img_diff, threshold=threshold, fail_fast=fail_fast)

        if mismatch == 0:
            return

        # Create new test_results folder
        failure_results_dir.mkdir(parents=True, exist_ok=True)
        img_diff.save(f"{failure_results_dir}/diff_{name}")
        img_a.save(f"{failure_results_dir}/actual_{name}")
        img_b.save(f"{failure_results_dir}/expected_{name}")

        # on ci, update the existing screenshots in place so we can download them
        if os.getenv("CI"):
            screenshot_file.write_bytes(img)

        # Still honor fail_fast if specifically requested
        if fail_fast:
            pytest.fail(f"{SNAPSHOT_MESSAGE_PREFIX} Snapshots DO NOT match! {name}")

        failures.append(f"{SNAPSHOT_MESSAGE_PREFIX} Snapshots DO NOT match! {name}")

    # Register finalizer to report all failures at the end of the test
    def finalize():
        if failures:
            pytest.fail("\n".join(failures))

    request.addfinalizer(finalize)

    return compare
