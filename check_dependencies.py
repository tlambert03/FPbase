"""Check release dates of project dependencies.

This script reads dependencies from pyproject.toml and queries PyPI to get:
- The release date of the currently specified version
- The release date of the most recent version (if different)
"""

from __future__ import annotations

import re
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

import requests


def parse_dependency(dep_spec: str) -> tuple[str, str | None]:
    """Parse a dependency specification into package name and version constraint.

    Parameters
    ----------
    dep_spec : str
        Dependency specification (e.g., 'django>=5.2,<5.3' or 'redis==5.0.4')

    Returns
    -------
    tuple[str, str | None]
        Package name and version constraint (or None if no constraint)
    """
    # Handle platform markers
    if ";" in dep_spec:
        dep_spec = dep_spec.split(";")[0].strip()

    # Extract package name (before extras or version specs)
    # Handle extras like package[extra1,extra2]
    match = re.match(r"([a-zA-Z0-9_\-]+)(\[.*?\])?(.*)", dep_spec)
    if not match:
        return dep_spec, None

    package_name = match.group(1).strip()
    # Extras are captured but not included in package name
    version_spec = match.group(3).strip() if match.group(3) else None

    return package_name, version_spec


def get_pinned_version(version_spec: str | None) -> str | None:
    """Extract pinned version from version specification.

    Parameters
    ----------
    version_spec : str | None
        Version specification (e.g., '>=5.2,<5.3' or '==5.0.4')

    Returns
    -------
    str | None
        Pinned version if using ==, otherwise None
    """
    if not version_spec:
        return None

    # Look for exact pin (==)
    match = re.search(r"==([0-9\.]+)", version_spec)
    return match.group(1) if match else None


def get_package_info(package_name: str) -> dict[str, Any] | None:
    """Fetch package information from PyPI.

    Parameters
    ----------
    package_name : str
        Name of the package

    Returns
    -------
    dict[str, Any] | None
        Package information from PyPI, or None if request fails
    """
    url = f"https://pypi.org/pypi/{package_name}/json"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"Error fetching {package_name}: {e}")
        return None


def format_date(date_str: str) -> str:
    """Format ISO date string to readable format.

    Parameters
    ----------
    date_str : str
        ISO format date string

    Returns
    -------
    str
        Formatted date string
    """
    dt = datetime.fromisoformat(date_str.replace("Z", "+00:00"))
    return dt.strftime("%Y-%m-%d")


def check_dependency(dep_spec: str) -> dict[str, Any]:
    """Check and gather release information for a dependency.

    Parameters
    ----------
    dep_spec : str
        Dependency specification from pyproject.toml

    Returns
    -------
    dict[str, Any]
        Dictionary containing package information and release dates
    """
    package_name, version_spec = parse_dependency(dep_spec)
    pinned_version = get_pinned_version(version_spec)

    result = {
        "package_name": package_name,
        "version_spec": version_spec,
        "pinned_version": pinned_version,
        "pinned_date": None,
        "latest_version": None,
        "latest_date": None,
        "error": None,
    }

    package_info = get_package_info(package_name)
    if not package_info:
        result["error"] = "Failed to fetch package info"
        return result

    latest_version = package_info["info"]["version"]
    releases = package_info.get("releases", {})

    result["latest_version"] = latest_version

    # Check pinned version
    if pinned_version and pinned_version in releases:
        release_info = releases[pinned_version]
        if release_info:
            result["pinned_date"] = format_date(release_info[0]["upload_time"])

    # Check latest version
    if latest_version in releases:
        latest_info = releases[latest_version]
        if latest_info:
            result["latest_date"] = format_date(latest_info[0]["upload_time"])

    return result


def print_result(result: dict[str, Any]) -> None:
    """Print formatted result for a dependency.

    Parameters
    ----------
    result : dict[str, Any]
        Result dictionary from check_dependency
    """
    print(f"\n{result['package_name']}")
    print(f"  Spec: {result['version_spec'] or 'any version'}")

    if result["error"]:
        print(f"  Error: {result['error']}")
        return

    if result["pinned_version"] and result["pinned_date"]:
        print(f"  Current version: {result['pinned_version']} (released {result['pinned_date']})")

    if result["latest_version"]:
        if result["latest_date"]:
            print(f"  Latest version:  {result['latest_version']} (released {result['latest_date']})")
        else:
            print(f"  Latest version:  {result['latest_version']}")

        if result["pinned_version"] and result["pinned_version"] != result["latest_version"]:
            print(f"  ⚠️  Update available: {result['pinned_version']} → {result['latest_version']}")


def main() -> None:
    """Main entry point for the script."""
    pyproject_path = Path("pyproject.toml")

    if not pyproject_path.exists():
        print("Error: pyproject.toml not found")
        return

    with open(pyproject_path, "rb") as f:
        pyproject = tomllib.load(f)

    dependencies = pyproject.get("project", {}).get("dependencies", [])

    if not dependencies:
        print("No dependencies found in pyproject.toml")
        return

    print(f"Checking {len(dependencies)} dependencies...\n")

    # Fetch all dependency information in parallel
    results = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_dep = {executor.submit(check_dependency, dep): dep for dep in dependencies}

        for future in as_completed(future_to_dep):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                dep = future_to_dep[future]
                print(f"Error processing {dep}: {e}")

    # Sort results by package name for consistent output
    results.sort(key=lambda r: r["package_name"].lower())

    # Print all results
    for result in results:
        print_result(result)

    # Gather packages with updates available and stale packages
    updates_available = []
    stale_packages = []

    for result in results:
        if result["error"]:
            continue

        # Check for updates
        if (
            result["pinned_version"]
            and result["latest_version"]
            and result["pinned_version"] != result["latest_version"]
        ):
            updates_available.append(
                {
                    "name": result["package_name"],
                    "current": result["pinned_version"],
                    "latest": result["latest_version"],
                    "current_date": result["pinned_date"],
                    "latest_date": result["latest_date"],
                }
            )

        # Check for stale packages
        if result["latest_date"]:
            latest_date = datetime.strptime(result["latest_date"], "%Y-%m-%d")
            days_old = (datetime.now() - latest_date).days

            if days_old >= 730:  # 2 years
                stale_packages.append(
                    {
                        "name": result["package_name"],
                        "version": result["latest_version"],
                        "date": result["latest_date"],
                        "days_old": days_old,
                    }
                )

    # Print updates available summary
    if updates_available:
        print("\n" + "=" * 70)
        print(f"\n📦 Updates available ({len(updates_available)}):\n")

        for pkg in updates_available:
            print(f"  {pkg['name']}: {pkg['current']} → {pkg['latest']}")
            if pkg["current_date"] and pkg["latest_date"]:
                print(f"    {pkg['current_date']} → {pkg['latest_date']}")

    # Print stale packages summary
    if stale_packages:
        print("\n" + "=" * 70)
        print(f"\n⚠️  Packages with no release in 2+ years ({len(stale_packages)}):\n")
        stale_packages.sort(key=lambda p: p["days_old"], reverse=True)

        for pkg in stale_packages:
            years = pkg["days_old"] / 365.25
            print(f"  {pkg['name']}: {pkg['version']} (last release {pkg['date']}, {years:.1f} years ago)")

    # Print summary if everything is up to date
    if not updates_available and not stale_packages:
        print("\n" + "=" * 70)
        print("\n✓ All packages are up to date and recently maintained")


if __name__ == "__main__":
    main()
