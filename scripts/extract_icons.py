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
import json
from pathlib import Path

from icon_extractors.fontawesome import FontAwesomeExtractor
from icon_extractors.lucide import LucideExtractor

# Map FPbase semantic icon names to library-specific names
# Format: {"fpbase-name": ("fa-style", "fa-icon-name")}
FONTAWESOME_ICON_MAP = {
    # UI & Navigation
    "info": ("fas", "info"),
    "warning": ("fas", "exclamation-circle"),
    "alert": ("fas", "exclamation-triangle"),
    "help": ("fas", "info-circle"),
    "question": ("fas", "question-circle"),
    "close": ("fas", "times"),
    "remove": ("fas", "times-circle"),
    "menu": ("fas", "list"),
    "grid": ("fas", "th"),
    "search": ("fas", "search"),
    "filter": ("fas", "filter"),
    "view": ("fas", "eye"),
    "settings": ("fas", "cog"),
    "edit": ("fas", "edit"),
    "delete": ("fas", "trash-alt"),
    "trash": ("fas", "trash"),
    "undo": ("fas", "undo"),
    "check": ("fas", "check"),
    "success": ("fas", "check-circle"),
    "selected": ("far", "check-square"),
    "unselected": ("far", "square"),
    "heart": ("fas", "heart"),
    # Actions
    "add": ("fas", "plus"),
    "add-item": ("fas", "plus-circle"),
    "download": ("fas", "download"),
    "upload": ("fas", "upload"),
    "share": ("fas", "share"),
    "share-square": ("fas", "share-square"),
    "link": ("fas", "link"),
    "external-link": ("fas", "external-link-alt"),
    "exchange": ("fas", "exchange-alt"),
    # Content
    "book": ("fas", "book"),
    "collection": ("fas", "book"),
    "quote": ("fas", "quote-left"),
    "photo": ("fas", "camera"),
    "chart": ("fas", "chart-area"),
    "table": ("fas", "table"),
    "flag": ("fas", "flag"),
    "flag-outline": ("far", "flag"),
    # Time & Status
    "clock": ("fas", "clock"),
    "spinner": ("fas", "spinner"),
    "lightbulb": ("fas", "lightbulb"),
    "sun": ("fas", "sun"),
    # Communication
    "email": ("fas", "envelope"),
    # Tools
    "wrench": ("fas", "wrench"),
    "keyboard": ("far", "keyboard"),
    # Social/External
    "google": ("fab", "google"),
    "twitter": ("fab", "x-twitter"),
    "orcid": ("fab", "orcid"),
}

# Map FPbase semantic icon names to Lucide icon names
LUCIDE_ICON_MAP = {
    # UI & Navigation
    "info": "info",
    "warning": "alert-circle",
    "alert": "alert-triangle",
    "help": "help-circle",
    "question": "help-circle",
    "close": "x",
    "remove": "x-circle",
    "menu": "list",
    "grid": "grid-3x3",
    "search": "search",
    "filter": "filter",
    "view": "eye",
    "settings": "settings",
    "edit": "edit",
    "delete": "trash-2",
    "trash": "trash",
    "undo": "undo",
    "check": "check",
    "success": "check-circle",
    "selected": "check-square",
    "unselected": "square",
    "heart": "heart",
    # Actions
    "add": "plus",
    "add-item": "plus-circle",
    "download": "download",
    "upload": "upload",
    "share": "share",
    "share-square": "share-2",
    "link": "link",
    "external-link": "external-link",
    "exchange": "arrow-left-right",
    # Content
    "book": "book",
    "collection": "book",
    "quote": "quote",
    "photo": "camera",
    "chart": "area-chart",
    "table": "table",
    "flag": "flag",
    "flag-outline": "flag",
    # Time & Status
    "clock": "clock",
    "spinner": "loader",
    "lightbulb": "lightbulb",
    "sun": "sun",
    # Communication
    "email": "mail",
    # Tools
    "wrench": "wrench",
    "keyboard": "keyboard",
    # Social/External (Lucide doesn't have brand icons, these will be missing)
    "google": None,
    "twitter": None,
    "orcid": None,
}


def extract_fontawesome_icons(output_dir: Path) -> None:
    """Extract FontAwesome icons to JSON."""
    print("Extracting FontAwesome icons...")
    extractor = FontAwesomeExtractor()

    icons = {}
    missing = []

    for fpbase_name, (style, fa_name) in FONTAWESOME_ICON_MAP.items():
        icon_data = extractor.extract_icon(fa_name, style)
        if icon_data:
            icons[fpbase_name] = icon_data
            print(f"  ✓ {fpbase_name} ({style} {fa_name})")
        else:
            missing.append(f"{fpbase_name} ({style} {fa_name})")
            print(f"  ✗ {fpbase_name} ({style} {fa_name}) - NOT FOUND")

    if missing:
        print(f"\nWarning: {len(missing)} icons not found:")
        for m in missing:
            print(f"  - {m}")

    # Save to JSON
    output_file = output_dir / "fontawesome.json"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as f:
        json.dump(icons, f, indent=2)

    print(f"\n✓ Extracted {len(icons)} icons to {output_file}")


def extract_lucide_icons(output_dir: Path) -> None:
    """Extract Lucide icons to JSON."""
    print("Extracting Lucide icons...")
    extractor = LucideExtractor()

    icons = {}
    missing = []

    for fpbase_name, lucide_name in LUCIDE_ICON_MAP.items():
        if lucide_name is None:
            print(f"  - {fpbase_name} - SKIPPED (no Lucide equivalent)")
            continue

        icon_data = extractor.extract_icon(lucide_name)
        if icon_data:
            icons[fpbase_name] = icon_data
            print(f"  ✓ {fpbase_name} ({lucide_name})")
        else:
            missing.append(f"{fpbase_name} ({lucide_name})")
            print(f"  ✗ {fpbase_name} ({lucide_name}) - NOT FOUND")

    if missing:
        print(f"\nWarning: {len(missing)} icons not found:")
        for m in missing:
            print(f"  - {m}")

    # Save to JSON
    output_file = output_dir / "lucide.json"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as f:
        json.dump(icons, f, indent=2)

    print(f"\n✓ Extracted {len(icons)} icons to {output_file}")


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
        extract_fontawesome_icons(args.output_dir)

    if args.library in ("lucide", "all"):
        extract_lucide_icons(args.output_dir)

    print("\n✓ Done!")


if __name__ == "__main__":
    main()
