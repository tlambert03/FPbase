"""A gunicorn worker that grew too large is replaced after its request."""

from __future__ import annotations

from types import SimpleNamespace
from typing import TYPE_CHECKING
from unittest.mock import Mock

from config import gunicorn

if TYPE_CHECKING:
    import pytest


def _post_request(rss_mb: float, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    monkeypatch.setattr(gunicorn, "peak_rss_mb", lambda: rss_mb)
    worker = SimpleNamespace(alive=True, pid=123, log=Mock())
    req = SimpleNamespace(method="POST", path="/graphql/")
    gunicorn.post_request(worker, req, {}, None)
    return worker


def test_large_worker_is_restarted(monkeypatch: pytest.MonkeyPatch) -> None:
    worker = _post_request(gunicorn.MAX_WORKER_RSS_MB + 1, monkeypatch)
    assert worker.alive is False
    worker.log.warning.assert_called_once()


def test_normal_worker_keeps_running(monkeypatch: pytest.MonkeyPatch) -> None:
    worker = _post_request(gunicorn.MAX_WORKER_RSS_MB - 1, monkeypatch)
    assert worker.alive is True
    worker.log.warning.assert_not_called()


def test_peak_rss_is_in_mb() -> None:
    # (a running Python process is tens of MB, not KB or GB)
    assert 10 < gunicorn.peak_rss_mb() < 10_000
