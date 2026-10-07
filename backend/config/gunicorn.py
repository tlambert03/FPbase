"""Gunicorn settings, loaded by the Procfile (which has the rest)."""

from __future__ import annotations

import os
import resource
import sys
from typing import Any

# A worker that grew past this (MB) is replaced after its request, rather than at
# --max-requests: Python keeps the memory of a large response (e.g. every dye's
# spectra) until the worker exits, and two such workers exceed the dyno's memory.
MAX_WORKER_RSS_MB = int(os.environ.get("GUNICORN_MAX_WORKER_RSS_MB", "250"))


def peak_rss_mb() -> float:
    """This process's peak resident memory, in MB."""
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # (bytes on macOS, KB on Linux)
    return maxrss / 1e6 if sys.platform == "darwin" else maxrss / 1e3


def post_request(worker: Any, req: Any, environ: dict, resp: Any) -> None:
    if (rss := peak_rss_mb()) > MAX_WORKER_RSS_MB:
        worker.log.warning(
            "Worker %s reached %d MB (over %d MB) by %s %s: restarting it",
            worker.pid,
            rss,
            MAX_WORKER_RSS_MB,
            req.method,
            req.path,
        )
        # (as --max-requests does: the worker exits after this request)
        worker.alive = False
