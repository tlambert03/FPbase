"""Lucide icon extractor.

Extracts icons from lucide-static npm package into the
standard FPbase icon format.
"""

from __future__ import annotations

from pathlib import Path
from xml.etree import ElementTree as ET

from . import IconData, IconExtractor, IconLibrary, IconPath

# Map FPbase semantic icon names to Lucide icon names
ICON_MAP = {
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
    "heart-outline": "heart",  # Lucide only has one heart (stroke-based)
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
    # Social/External (Lucide doesn't have brand icons)
    "google": None,
    "twitter": None,
    "orcid": None,
}

# Icons that should use filled variant
# Lucide's recommended workaround: https://lucide.dev/guide/advanced/filled-icons
# Uses fill="currentColor" and strokeWidth="0" instead of stroke
# Works best for icons with closed paths (heart, checkboxes, etc.)
FILLED_ICONS = {"heart", "selected"}

# Lucide default SVG attributes
DEFAULTS = {
    "stroke": "currentColor",
    "strokeWidth": "2",
    "strokeLinecap": "round",
    "strokeLinejoin": "round",
    "fill": "none",
}

# Lucide icons render slightly smaller than FontAwesome at 1em
# Increase scale slightly for better visual consistency
SCALE = "1.1em"


def circle_to_path(cx: float, cy: float, r: float) -> str:
    """Convert SVG circle to path data.

    Parameters
    ----------
    cx, cy : float
        Circle center
    r : float
        Circle radius

    Returns
    -------
    str
        SVG path data
    """
    # Circle using two arc commands
    return f"M {cx - r} {cy} a {r} {r} 0 1 0 {r * 2} 0 a {r} {r} 0 1 0 {-r * 2} 0"


def line_to_path(x1: float, y1: float, x2: float, y2: float) -> str:
    """Convert SVG line to path data.

    Parameters
    ----------
    x1, y1 : float
        Start point
    x2, y2 : float
        End point

    Returns
    -------
    str
        SVG path data
    """
    return f"M {x1} {y1} L {x2} {y2}"


def rect_to_path(x: float, y: float, width: float, height: float, rx: float = 0, ry: float = 0) -> str:
    """Convert SVG rect to path data.

    Parameters
    ----------
    x, y : float
        Rectangle position
    width, height : float
        Rectangle dimensions
    rx, ry : float
        Corner radius (for rounded rectangles)

    Returns
    -------
    str
        SVG path data
    """
    if rx == 0 and ry == 0:
        # Simple rectangle (no rounded corners)
        return f"M {x} {y} h {width} v {height} h {-width} Z"

    # Rounded rectangle
    # Clamp radius to half the smallest dimension
    rx = min(rx, width / 2)
    ry = min(ry, height / 2) if ry > 0 else rx

    return (
        f"M {x + rx} {y} "
        f"h {width - 2 * rx} "
        f"a {rx} {ry} 0 0 1 {rx} {ry} "
        f"v {height - 2 * ry} "
        f"a {rx} {ry} 0 0 1 {-rx} {ry} "
        f"h {-(width - 2 * rx)} "
        f"a {rx} {ry} 0 0 1 {-rx} {-ry} "
        f"v {-(height - 2 * ry)} "
        f"a {rx} {ry} 0 0 1 {rx} {-ry} "
        f"Z"
    )


class LucideExtractor(IconExtractor):
    """Extract icons from Lucide.

    Parameters
    ----------
    base_path : Path | None
        Path to lucide-static package.
        If None, uses default node_modules location.
    """

    def __init__(self, base_path: Path | None = None):
        if base_path is None:
            # Default to node_modules location
            base_path = Path(__file__).parent.parent.parent / "node_modules" / "lucide-static"

        self.base_path = base_path
        self.icons_path = base_path / "icons"

        if not self.icons_path.exists():
            raise FileNotFoundError(
                f"Lucide icons not found at {self.icons_path}. Make sure lucide-static is installed."
            )

    def extract_icon(self, icon_name: str, filled: bool = False) -> IconData | None:
        """Extract a Lucide icon.

        Parameters
        ----------
        icon_name : str
            The Lucide icon name (e.g., "external-link", "heart")
        filled : bool
            If True, modify the icon to be filled (set fill="currentColor", stroke-width="0")

        Returns
        -------
        IconData | None
            The icon data in standard format, or None if not found
        """
        svg_path = self.icons_path / f"{icon_name}.svg"

        if not svg_path.exists():
            return None

        try:
            # Parse the SVG file
            tree = ET.parse(svg_path)
            root = tree.getroot()

            # Extract viewBox
            viewBox = root.get("viewBox", "0 0 24 24")

            # Extract paths (just the path data)
            paths: list[IconPath] = []

            # Process <path> elements
            for path_elem in root.findall(".//{http://www.w3.org/2000/svg}path"):
                path_data: IconPath = {
                    "d": path_elem.get("d", ""),
                }

                # Handle filled variant (Lucide workaround from their docs)
                if filled:
                    # For filled variant: override with fill at path level
                    path_data["fill"] = "currentColor"  # type: ignore
                    path_data["strokeWidth"] = "0"  # type: ignore

                # Remove None values
                path_data = {k: v for k, v in path_data.items() if v is not None}  # type: ignore
                paths.append(path_data)  # type: ignore

            # Process <rect> elements (convert to path)
            for rect_elem in root.findall(".//{http://www.w3.org/2000/svg}rect"):
                x = float(rect_elem.get("x", "0"))
                y = float(rect_elem.get("y", "0"))
                width = float(rect_elem.get("width", "0"))
                height = float(rect_elem.get("height", "0"))
                rx = float(rect_elem.get("rx", "0"))
                ry = float(rect_elem.get("ry", rx))  # ry defaults to rx if not specified

                path_d = rect_to_path(x, y, width, height, rx, ry)
                path_data: IconPath = {"d": path_d}

                # Apply same fill logic as for path elements
                if filled:
                    path_data["fill"] = "currentColor"  # type: ignore
                    path_data["strokeWidth"] = "0"  # type: ignore
                else:
                    rect_fill = rect_elem.get("fill")
                    if rect_fill and rect_fill != "none":
                        path_data["fill"] = rect_fill  # type: ignore

                paths.append(path_data)  # type: ignore

            # Process <circle> elements (convert to path)
            for circle_elem in root.findall(".//{http://www.w3.org/2000/svg}circle"):
                cx = float(circle_elem.get("cx", "0"))
                cy = float(circle_elem.get("cy", "0"))
                r = float(circle_elem.get("r", "0"))

                path_d = circle_to_path(cx, cy, r)
                path_data: IconPath = {"d": path_d}

                # Apply same fill logic as for path elements
                if filled:
                    path_data["fill"] = "currentColor"  # type: ignore
                    path_data["strokeWidth"] = "0"  # type: ignore
                else:
                    circle_fill = circle_elem.get("fill")
                    if circle_fill and circle_fill != "none":
                        path_data["fill"] = circle_fill  # type: ignore

                paths.append(path_data)  # type: ignore

            # Process <line> elements (convert to path)
            for line_elem in root.findall(".//{http://www.w3.org/2000/svg}line"):
                x1 = float(line_elem.get("x1", "0"))
                y1 = float(line_elem.get("y1", "0"))
                x2 = float(line_elem.get("x2", "0"))
                y2 = float(line_elem.get("y2", "0"))

                path_d = line_to_path(x1, y1, x2, y2)
                path_data: IconPath = {"d": path_d}

                # Apply same fill logic as for path elements
                if filled:
                    path_data["fill"] = "currentColor"  # type: ignore
                    path_data["strokeWidth"] = "0"  # type: ignore
                else:
                    line_fill = line_elem.get("fill")
                    if line_fill and line_fill != "none":
                        path_data["fill"] = line_fill  # type: ignore

                paths.append(path_data)  # type: ignore

            # Build icon data
            icon_data = IconData(viewBox=viewBox, paths=paths)

            if filled:
                # For filled variant, only store the overrides (fill and strokeWidth)
                # These differ from Lucide defaults and need to be explicit
                icon_data["fill"] = "currentColor"  # type: ignore
                icon_data["strokeWidth"] = "0"  # type: ignore
            else:
                # For outline variant, only store attributes that differ from defaults
                # Lucide defaults: stroke="currentColor", strokeWidth="2",
                #                  strokeLinecap="round", strokeLinejoin="round", fill="none"
                stroke = root.get("stroke", "currentColor")
                if stroke != "currentColor":
                    icon_data["stroke"] = stroke  # type: ignore

                stroke_width = root.get("stroke-width", "2")
                if stroke_width != "2":
                    icon_data["strokeWidth"] = stroke_width  # type: ignore

                stroke_linecap = root.get("stroke-linecap", "round")
                if stroke_linecap != "round":
                    icon_data["strokeLinecap"] = stroke_linecap  # type: ignore

                stroke_linejoin = root.get("stroke-linejoin", "round")
                if stroke_linejoin != "round":
                    icon_data["strokeLinejoin"] = stroke_linejoin  # type: ignore

                # Don't store fill="none" since it's the default

            return icon_data

        except Exception as e:
            print(f"Error parsing {svg_path}: {e}")
            return None

    def get_available_icons(self) -> list[str]:
        """Get all available Lucide icon names.

        Returns
        -------
        list[str]
            List of available icon names
        """
        if not self.icons_path.exists():
            return []

        # Get all .svg files and strip the extension
        return [f.stem for f in self.icons_path.glob("*.svg")]

    def extract_library(self) -> IconLibrary:
        """Extract all Lucide icons defined in ICON_MAP.

        Returns
        -------
        IconLibrary
            Complete library data ready to be written to JSON
        """
        print("Extracting Lucide icons...")

        icons = {}
        missing = []

        for fpbase_name, lucide_name in ICON_MAP.items():
            if lucide_name is None:
                print(f"  - {fpbase_name} - SKIPPED (no Lucide equivalent)")
                continue

            # Use filled variant for specific icons
            filled = fpbase_name in FILLED_ICONS
            icon_data = self.extract_icon(lucide_name, filled=filled)
            if icon_data:
                icons[fpbase_name] = icon_data
                variant = " [filled]" if filled else ""
                print(f"  ✓ {fpbase_name} ({lucide_name}){variant}")
            else:
                missing.append(f"{fpbase_name} ({lucide_name})")
                print(f"  ✗ {fpbase_name} ({lucide_name}) - NOT FOUND")

        if missing:
            print(f"\nWarning: {len(missing)} icons not found:")
            for m in missing:
                print(f"  - {m}")

        return IconLibrary(icons=icons, scale=SCALE, defaults=DEFAULTS)
