"""FontAwesome icon extractor.

Extracts icons from @fortawesome/fontawesome-free npm package into the
standard FPbase icon format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal
from xml.etree import ElementTree as ET

from . import IconData, IconExtractor, IconPath


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
                    "fill": path_elem.get("fill"),
                    "fillRule": path_elem.get("fill-rule"),
                    "clipRule": path_elem.get("clip-rule"),
                }
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
