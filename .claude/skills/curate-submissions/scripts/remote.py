#!/usr/bin/env python3
"""Run a script from this folder inside Django on a Heroku one-off dyno.

The script is piped over stdin to `heroku run ... python -`, so nothing needs to be
deployed.  `apply` is a dry run (rolled back) unless `--commit` is passed.
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from views import annotate, load_views

HERE = Path(__file__).parent
APP = os.environ.get("FPBASE_HEROKU_APP", "fpbase")
TIMEOUT = int(os.environ.get("FPBASE_REMOTE_TIMEOUT", "300"))


def _text(x: str | bytes | None) -> str:
    return x.decode(errors="replace") if isinstance(x, bytes) else (x or "")


def run_remote(script: str, params: dict, *, write: bool = False) -> dict:
    b64 = base64.b64encode(json.dumps(params).encode()).decode()
    source = "\n".join(
        [
            (HERE / "_bootstrap.py").read_text(),
            # everything except `apply` runs in a session where Postgres refuses writes
            "" if write else "read_only()",
            f"import base64; PARAMS = json.loads(base64.b64decode('{b64}'))",
            (HERE / script).read_text(),
        ]
    )
    # --no-notify: otherwise the CLI pops a desktop "dyno is up" notification on every call
    cmd = ["heroku", "run", "--no-tty", "--no-notify", "--exit-code", "-a", APP, "--"]
    cmd += ["python", "-"]
    try:
        proc = subprocess.run(cmd, input=source, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired as e:
        # a hung `apply` would keep a transaction (and row locks) open on production:
        # stop our own one-off dyno rather than just abandoning it
        out, err = (_text(e.stdout), _text(e.stderr))
        for dyno in set(re.findall(r"run\.\d+", err)):
            subprocess.run(["heroku", "ps:stop", dyno, "-a", APP], capture_output=True)
            print(f"timed out after {TIMEOUT}s: stopped {dyno}", file=sys.stderr)
        sys.exit(f"remote script timed out; nothing was committed.\n{out}\n{err}")
    out = proc.stdout
    if "<<<JSON" not in out or "JSON>>>" not in out:
        sys.exit(f"remote script failed (exit {proc.returncode}):\n{out}\n{proc.stderr}")
    return json.loads(out.split("<<<JSON")[1].split("JSON>>>")[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    fetch = sub.add_parser("fetch", help="dump pending submissions (read-only)")
    fetch.add_argument("--kind", choices=["proteins", "spectra", "all"], default="all")
    fetch.add_argument("--limit", type=int)
    fetch.add_argument("--offset", type=int, default=0)
    fetch.add_argument("--slugs", nargs="+", help="only these protein slugs")
    fetch.add_argument("--summary", action="store_true", help="one line per protein")
    fetch.add_argument("-o", "--output", type=Path)

    triage = sub.add_parser("triage", help="net change of every pending protein (read-only)")
    triage.add_argument("--slugs", nargs="+")
    triage.add_argument("-o", "--output", type=Path)

    audit = sub.add_parser("audit", help="machine-checkable health of proteins (read-only)")
    audit.add_argument("--top", type=int, help="the N most-viewed protein pages")
    audit.add_argument("--slugs", nargs="+")
    audit.add_argument("-o", "--output", type=Path)

    prov = sub.add_parser("provenance", help="who set each value, and when (read-only)")
    prov.add_argument("slug")
    prov.add_argument("--around", type=int, nargs="+", help="print spectrum data around these nm")
    prov.add_argument("-o", "--output", type=Path)

    apply = sub.add_parser("apply", help="apply a decisions file (dry run by default)")
    apply.add_argument("decisions", type=Path)
    apply.add_argument("--commit", action="store_true", help="actually write to production")
    apply.add_argument("--moderator", default=os.environ.get("FPBASE_MODERATOR"))
    apply.add_argument("-o", "--output", type=Path)

    args = parser.parse_args()
    if args.cmd == "triage":
        result = run_remote("triage_pending.py", {"slugs": args.slugs})
        # most-viewed first (no-op until `views.py ingest` has been run)
        result = annotate(result, load_views())
    elif args.cmd == "audit":
        slugs = args.slugs or list(load_views())[: args.top or 100]
        result = annotate(run_remote("audit_proteins.py", {"slugs": slugs}), load_views())
    elif args.cmd == "provenance":
        result = run_remote("provenance.py", {"slug": args.slug, "around": args.around})
    elif args.cmd == "fetch":
        params = {k: getattr(args, k) for k in ("kind", "limit", "offset", "slugs", "summary")}
        result = run_remote("fetch_pending.py", params)
    else:
        if not args.moderator:
            sys.exit("pass --moderator <fpbase staff username> (or set FPBASE_MODERATOR)")
        params = {
            "decisions": json.loads(args.decisions.read_text()),
            "commit": args.commit,
            "moderator": args.moderator,
        }
        result = run_remote("apply_decisions.py", params, write=True)

    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text)
        print(f"wrote {args.output}")
    else:
        print(text)


if __name__ == "__main__":
    main()
