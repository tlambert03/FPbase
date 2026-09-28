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


def test_request_finished_includes_duration() -> None:
    request = RequestFactory().get("/")
    signals.bind_extra_request_metadata.send(
        sender=None, request=request, logger=None, log_kwargs={}
    )
    log_kwargs: dict = {}
    signals.bind_extra_request_finished_metadata.send(
        sender=None, request=request, logger=None, response=None, log_kwargs=log_kwargs
    )
    assert isinstance(log_kwargs["duration_ms"], int)


@pytest.mark.django_db
def test_graphql_bad_request_is_logged(client, caplog) -> None:
    query = "{ protein(id: 1) { notAField } }"
    with caplog.at_level("WARNING", logger="fpbase.views"):
        response = client.post("/graphql/", {"query": query}, content_type="application/json")
    assert response.status_code == 400
    (record,) = [r for r in caplog.records if r.getMessage() == "GraphQL bad request"]
    assert record.query == query
    assert any("notAField" in e for e in record.errors)
