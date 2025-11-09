"""FontAwesome icon extractor.

Extracts icons from @fortawesome/fontawesome-free npm package into the
standard FPbase icon format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal
from xml.etree import ElementTree as ET

from . import IconData, IconExtractor, IconLibrary, IconPath

# Map FPbase semantic icon names to FontAwesome-specific names
# Format: {"fpbase-name": ("fa-style", "fa-icon-name")}
ICON_MAP = {
    # UI & Navigation
    "info": ("fas", "info"),
    "info-circle": ("fas", "info-circle"),
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
    "heart-outline": ("far", "heart"),
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

# FontAwesome uses fill-based icons, so no SVG-level defaults needed
DEFAULTS = {"fill": "currentColor"}


class FontAwesomeExtractor(IconExtractor):
    """Extract icons from FontAwesome Free.

    Parameters
    ----------
    base_path : Path | None
        Path to @fortawesome/fontawesome-free package.
        If None, uses default node_modules location.
    """

    STYLE_MAP = {
        "fas": "solid",
        "far": "regular",
        "fab": "brands",
    }

    def __init__(self, base_path: Path | None = None):
        if base_path is None:
            # Default to node_modules location
            base_path = Path(__file__).parent.parent.parent / "node_modules" / "@fortawesome" / "fontawesome-free"

        self.base_path = base_path
        self.svgs_path = base_path / "svgs"

        if not self.svgs_path.exists():
            raise FileNotFoundError(
                f"FontAwesome SVGs not found at {self.svgs_path}. "
                "Make sure @fortawesome/fontawesome-free is installed."
            )

    def extract_icon(self, icon_name: str, style: Literal["fas", "far", "fab"] = "fas") -> IconData | None:
        """Extract a FontAwesome icon.

        Parameters
        ----------
        icon_name : str
            The FontAwesome icon name (e.g., "external-link-alt", "heart")
        style : Literal["fas", "far", "fab"]
            The FontAwesome style: fas (solid), far (regular), fab (brands)

        Returns
        -------
        IconData | None
            The icon data in standard format, or None if not found
        """
        style_dir = self.STYLE_MAP.get(style)
        if not style_dir:
            print(f"Warning: Unknown FontAwesome style '{style}'")
            return None

        svg_path = self.svgs_path / style_dir / f"{icon_name}.svg"

        if not svg_path.exists():
            return None

        try:
            # Parse the SVG file
            tree = ET.parse(svg_path)
            root = tree.getroot()

            # Extract viewBox
            viewBox = root.get("viewBox", "0 0 512 512")

            # Extract paths
            paths: list[IconPath] = []
            for path_elem in root.findall(".//{http://www.w3.org/2000/svg}path"):
                path_data: IconPath = {
                    "d": path_elem.get("d", ""),
                }
                for attr in ["fill", "fillRule", "clipRule"]:
                    val = path_elem.get(attr.replace("R", "-r").lower())
                    if val is not None and val != DEFAULTS.get(attr):
                        path_data[attr] = val

                # Remove None values
                path_data = {k: v for k, v in path_data.items() if v is not None}  # type: ignore
                paths.append(path_data)  # type: ignore

            return IconData(viewBox=viewBox, paths=paths)

        except Exception as e:
            print(f"Error parsing {svg_path}: {e}")
            return None

    def get_available_icons(self, style: Literal["fas", "far", "fab"] = "fas") -> list[str]:
        """Get all available icon names for a given style.

        Parameters
        ----------
        style : Literal["fas", "far", "fab"]
            The FontAwesome style

        Returns
        -------
        list[str]
            List of available icon names
        """
        style_dir = self.STYLE_MAP.get(style)
        if not style_dir:
            return []

        style_path = self.svgs_path / style_dir
        if not style_path.exists():
            return []

        # Get all .svg files and strip the extension
        return [f.stem for f in style_path.glob("*.svg")]

    def extract_library(self) -> IconLibrary:
        """Extract all FontAwesome icons defined in ICON_MAP.

        Returns
        -------
        IconLibrary
            Complete library data ready to be written to JSON
        """
        print("Extracting FontAwesome icons...")

        icons = {}
        missing = []

        for fpbase_name, (style, fa_name) in ICON_MAP.items():
            icon_data = self.extract_icon(fa_name, style)
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

        return IconLibrary(icons=icons, defaults=DEFAULTS)
