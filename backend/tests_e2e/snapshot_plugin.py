"""Pytest plugin for visual snapshot testing with Playwright.

This plugin provides fixtures and configuration for visual regression testing.
It's registered in conftest.py via pytest_plugins to make fixtures auto-discoverable.
"""

from __future__ import annotations

import math
import os
import shutil
import sys
from contextlib import suppress
from copy import copy
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

if TYPE_CHECKING:
    from collections.abc import Callable

DEFAULT_MASK_SELECTORS = [".highcharts-legend", "img[src*='GFP_spinner']"]


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add custom command line options for visual snapshot tests."""
    parser.addoption(
        "--assert-snapshots",
        action="store_true",
        default=False,
        help="Enable visual snapshot testing (disabled by default, never runs on CI)",
    )
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Update visual snapshots instead of comparing (for use with --assert-snapshots)",
    )
    parser.addoption(
        "--snapshot-threshold",
        type=float,
        default=0.1,
        help="Threshold for pixel differences in visual snapshot comparisons (default: 0.1)",
    )
    parser.addoption(
        "--min-percent-diff",
        type=float,
        default=0.1,  # 0.1 %
        help="Minimum percent pixels allowed to be different before "
        "failing a visual snapshot test (default: 0) out of 100",
    )


def _visual_snapshots_enabled(config: pytest.Config) -> bool:
    """Check if visual snapshots are enabled via CLI flag or environment variable."""
    return (
        config.getoption("--assert-snapshots", False)
        or config.getoption("--update-snapshots", False)
        or os.environ.get("VISUAL_SNAPSHOTS", "").lower()
        in (
            "1",
            "true",
            "yes",
        )
    )


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

    root_dir = Path(pytestconfig.rootpath)  # type: ignore[attr-defined]

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

        noop.NOOP = True  # type: ignore[attr-defined]
        return noop

    from io import BytesIO

    import pytest
    from PIL import Image
    from pixelmatch.contrib.PIL import pixelmatch

    test_function_name = request.node.name
    SNAPSHOT_MESSAGE_PREFIX = "[playwright-assert-snapshot]"

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
    global_snapshot_threshold = pytestconfig.getoption("--snapshot-threshold", 0.1)
    min_percent_diff = pytestconfig.getoption("--min-percent-diff", 0)

    mask_selectors: list[str] = copy(DEFAULT_MASK_SELECTORS)
    update_snapshot = pytestconfig.getoption("--update-snapshots", False)

    # for automatically naming multiple assertions
    counter = 0
    # Collection to store failures
    failures = []

    def _create_locators_from_selectors(page, selectors: list[str]):
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
        if hasattr(img_or_page, "screenshot"):  # Duck-type check for Page
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

        # Calculate percentage
        mismatch_percentage = (mismatch / math.prod(img_a.size)) * 100
        if mismatch_percentage <= min_percent_diff:
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
        msg = (
            f"{SNAPSHOT_MESSAGE_PREFIX} Snapshots DO NOT match! {name}. "
            f"Mismatched pixels: {mismatch} ({mismatch_percentage:.2f}%)"
        )
        if fail_fast:
            pytest.fail(msg)

        failures.append(msg)

    # Register finalizer to report all failures at the end of the test
    def finalize():
        if failures:
            first_line = f"{len(failures)} visual snapshot test(s) failed:\n"
            pytest.fail(first_line + "\n".join(failures))

    if pytestconfig.getoption("--assert-snapshots", False):
        request.addfinalizer(finalize)

    return compare
