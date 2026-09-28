"""Extra fields on django-structlog's per-request log lines."""

from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from django.dispatch import receiver
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
