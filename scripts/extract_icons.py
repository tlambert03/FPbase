#!/usr/bin/env python3
"""Extract icons from various libraries into FPbase's standard format.

This script maps FPbase's semantic icon names to library-specific icon names
and extracts the SVG data into JSON files that can be used by the Django
template tags.

Usage:
    python scripts/extract_icons.py fontawesome
    python scripts/extract_icons.py lucide
    python scripts/extract_icons.py all
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from icon_extractors import fontawesome, lucide


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract icons from libraries into FPbase standard format")
    parser.add_argument(
        "library",
        choices=["fontawesome", "lucide", "all"],
        help="Icon library to extract (or 'all' for all libraries)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent.parent / "backend" / "fpbase" / "static" / "icons",
        help="Output directory for JSON files",
    )

    args = parser.parse_args()

    if args.library in ("fontawesome", "all"):
        fontawesome.FontAwesomeExtractor.extract_to(args.output_dir / "fontawesome.json")

    if args.library in ("lucide", "all"):
        lucide.LucideExtractor.extract_to(args.output_dir / "lucide.json")

    print("\n✓ Done!")


if __name__ == "__main__":
    main()
    subprocess.run(["prek", "run", "end-of-file-fixer", "-a"])
