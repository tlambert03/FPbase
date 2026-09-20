#!/usr/bin/env python3
"""Page views per protein, for deciding what deserves attention first.

`ingest` turns a saved Google Analytics `run_report` result (dimension `pagePath`, metric
`screenPageViews`) into `.curation/analytics/protein_views.json`.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path

VIEWS_FILE = Path(".curation/analytics/protein_views.json")
DETAIL_PAGE = re.compile(r"^/protein/([^/]+)/?$")  # not /history/, /bleach/, ...


def parse_report(report: dict) -> dict[str, int]:
    views: Counter[str] = Counter()
    for row in report["rows"]:
        if m := DETAIL_PAGE.match(row["dimension_values"][0]["value"]):
            views[m.group(1).lower()] += int(row["metric_values"][0]["value"])
    return dict(views.most_common())


def load_views(path: Path = VIEWS_FILE) -> dict[str, int]:
    return json.loads(path.read_text())["views"] if path.exists() else {}


def annotate(rows: list[dict], views: dict[str, int]) -> list[dict]:
    """Add `views_365d` and site-wide `rank` to each row; most viewed first."""
    rank = {slug: i for i, slug in enumerate(views, 1)}
    for row in rows:
        row["views_365d"] = views.get(row["slug"], 0)
        row["rank"] = rank.get(row["slug"])
    return sorted(rows, key=lambda r: -r["views_365d"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    ingest = sub.add_parser("ingest", help="GA run_report JSON -> protein_views.json")
    ingest.add_argument("report", type=Path)
    ingest.add_argument("-o", "--output", type=Path, default=VIEWS_FILE)
    args = parser.parse_args()

    views = parse_report(json.loads(args.report.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"fetched": str(date.today()), "views": views}, indent=0))
    total = sum(views.values())
    share = sum(list(views.values())[:100]) / total
    print(f"wrote {args.output}: {len(views)} pages, {total:,} views, top 100 = {share:.0%}")


if __name__ == "__main__":
    main()
