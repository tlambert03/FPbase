from __future__ import annotations

import pytest
from django.test import RequestFactory
from django_structlog import signals


@pytest.mark.parametrize(
    "signal",
    [
        signals.bind_extra_request_finished_metadata,
        signals.bind_extra_request_failed_metadata,
    ],
)
def test_request_log_lines_include_client_headers(signal) -> None:
    request = RequestFactory().get(
        "/", headers={"user-agent": "bot/1.0", "referer": "https://example.com/"}
    )
    log_kwargs = {"code": 200, "request": "GET /"}
    signal.send(
        sender=None,
        request=request,
        logger=None,
        response=None,
        exception=None,
        log_kwargs=log_kwargs,
    )
    assert log_kwargs["user_agent"] == "bot/1.0"
    assert log_kwargs["referer"] == "https://example.com/"
