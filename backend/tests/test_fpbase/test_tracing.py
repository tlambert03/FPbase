from __future__ import annotations

import pytest

from fpbase.tracing import traces_sampler

CHROME = "Mozilla/5.0 (Macintosh) AppleWebKit/537.36 Chrome/154.0.0.0 Safari/537.36"


def _ctx(path: str = "/protein/egfp/", ua: str = CHROME, parent: bool | None = None) -> dict:
    return {
        "wsgi_environ": {"PATH_INFO": path, "HTTP_USER_AGENT": ua},
        "parent_sampled": parent,
    }


@pytest.mark.parametrize(
    "ua",
    [
        "Mozilla/5.0 (compatible; Googlebot/2.1; +http://www.google.com/bot.html)",
        "Mozilla/5.0 AppleWebKit/537.36 (compatible; bingbot/2.0) Chrome/116",
        "Claude-User (claude-code/2.1.226; +https://support.anthropic.com/)",
        "GitHubCopilotRuntime-WebFetch",
        "Mozilla/5.0 HeadlessChrome/120.0.0.0",
    ],
)
def test_bots_not_traced(ua: str) -> None:
    assert traces_sampler(_ctx(ua=ua, parent=True)) == 0


def test_static_not_traced() -> None:
    assert traces_sampler(_ctx(path="/static/main.js", parent=True)) == 0


@pytest.mark.parametrize("ua", [CHROME, "python-requests/2.34.2", "fpbase-py/0.2.1"])
@pytest.mark.parametrize("parent", [True, False])
def test_follows_parent_decision(ua: str, parent: bool) -> None:
    assert traces_sampler(_ctx(ua=ua, parent=parent)) == float(parent)
