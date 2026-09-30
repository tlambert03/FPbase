"""Extra fields on django-structlog's per-request log lines."""

from __future__ import annotations

import json
import re
import time
from typing import TYPE_CHECKING, Any

from django.dispatch import receiver
from django.http.request import RawPostDataException
from django_structlog import signals

if TYPE_CHECKING:
    from django.http import HttpRequest


@receiver(signals.bind_extra_request_metadata)
def mark_request_start(request: HttpRequest, **kwargs: Any) -> None:
    request._log_start = time.monotonic()  # pyright: ignore[reportAttributeAccessIssue]


@receiver(signals.bind_extra_request_finished_metadata)
@receiver(signals.bind_extra_request_failed_metadata)
def add_client_headers(request: HttpRequest, log_kwargs: dict[str, Any], **kwargs: Any) -> None:
    # production only logs `request_finished`, so it carries what Heroku's router lacks
    log_kwargs["user_agent"] = request.headers.get("user-agent")
    log_kwargs["referer"] = request.headers.get("referer")
    # time spent in Django (from the structlog middleware on), so the line stands alone
    if (start := getattr(request, "_log_start", None)) is not None:
        log_kwargs["duration_ms"] = round((time.monotonic() - start) * 1000)
    if request.path.startswith("/graphql"):
        log_kwargs["graphql_operation"] = _graphql_operation(request)


# `query getMicroscope(...) {` -> getMicroscope; an anonymous `{ spectra {...} }` -> spectra
_OPERATION_RE = re.compile(
    r"^\s*(?:query|mutation|subscription)\s+(\w+)|^\s*(?:query\s*)?\{\s*(\w+)"
)


def _graphql_operation(request: HttpRequest) -> str | None:
    """Which GraphQL operation a request ran, so slow ones can be told apart."""
    try:
        payload = json.loads(request.body) if request.method == "POST" else request.GET
    except (ValueError, RawPostDataException):
        return None
    if not isinstance(payload, dict):  # a batch
        return None
    if name := payload.get("operationName"):
        return str(name)[:100]
    match = _OPERATION_RE.match(str(payload.get("query") or ""))
    return (match.group(1) or match.group(2)) if match else None
