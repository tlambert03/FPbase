"""Utility to get git version information for deployment tracking."""

import subprocess
from functools import lru_cache


@lru_cache(maxsize=1)
def get_git_sha() -> str:
    """
    Get the current git SHA.

    Returns the full 40-character SHA if available, otherwise a default message.
    The result is cached to avoid repeated git calls.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return "unknown"


@lru_cache(maxsize=1)
def get_git_sha_short() -> str:
    """
    Get the short version of the git SHA (first 7 characters).

    Returns the short SHA if available, otherwise a default message.
    The result is cached to avoid repeated git calls.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return "unknown"


@lru_cache(maxsize=1)
def get_github_url(sha: str | None = None) -> str:
    """
    Get the GitHub URL for the current commit.

    Args:
        sha: Optional SHA to use. If not provided, will get the current HEAD SHA.

    Returns:
        The GitHub URL for the commit, or an empty string if unable to determine.
    """
    if sha is None:
        sha = get_git_sha()

    if sha == "unknown":
        return ""

    # FPbase repository URL
    # TODO: Make this configurable if needed
    base_url = "https://github.com/tlambert03/FPbase"
    return f"{base_url}/commit/{sha}"
