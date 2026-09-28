"""Extra fields on django-structlog's per-request log lines."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from django.dispatch import receiver
from django_structlog import signals

if TYPE_CHECKING:
    from django.http import HttpRequest


@receiver(signals.bind_extra_request_finished_metadata)
@receiver(signals.bind_extra_request_failed_metadata)
def add_client_headers(request: HttpRequest, log_kwargs: dict[str, Any], **kwargs: Any) -> None:
    # production only logs `request_finished`, so it carries what Heroku's router lacks
    log_kwargs["user_agent"] = request.headers.get("user-agent")
    log_kwargs["referer"] = request.headers.get("referer")
