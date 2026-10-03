"""Which requests Sentry traces (`traces_sampler` in production settings)."""

from __future__ import annotations

import re
from typing import Any

import sentry_sdk

# crawlers and AI agents that name themselves; scripted API clients (python-requests,
# curl, fpbase-py) are real API usage and stay in. Keep in sync with sentry-init.js
BOT_UA = re.compile(
    r"bot|crawl|spider|slurp|scrap|headless|phantom|preview|facebookexternalhit"
    r"|claude|chatgpt|copilot|perplexity",
    re.IGNORECASE,
)
SKIP_PATHS = ("/static/", "/media/", "/robots.txt", "/favicon")


def traces_sampler(context: dict[str, Any]) -> float:
    if environ := context.get("wsgi_environ"):
        if environ.get("PATH_INFO", "").startswith(SKIP_PATHS):
            return 0
        if BOT_UA.search(environ.get("HTTP_USER_AGENT", "")):
            return 0
    # keep the browser's decision for its API calls, so its traces stay whole
    if (parent := context.get("parent_sampled")) is not None:
        return float(parent)
    return sentry_sdk.get_client().options["traces_sample_rate"] or 0
