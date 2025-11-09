"""Lucide icon extractor.

Extracts icons from lucide-static npm package into the
standard FPbase icon format.
"""

from __future__ import annotations

from pathlib import Path
from xml.etree import ElementTree as ET

from . import IconData, IconExtractor, IconPath


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

    def extract_icon(self, icon_name: str) -> IconData | None:
        """Extract a Lucide icon.

        Parameters
        ----------
        icon_name : str
            The Lucide icon name (e.g., "external-link", "heart")

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

            # Lucide uses stroke-based icons, so we need to extract stroke attributes
            stroke = root.get("stroke", "currentColor")
            stroke_width = root.get("stroke-width", "2")
            stroke_linecap = root.get("stroke-linecap", "round")
            stroke_linejoin = root.get("stroke-linejoin", "round")

            # Extract paths
            paths: list[IconPath] = []
            for path_elem in root.findall(".//{http://www.w3.org/2000/svg}path"):
                path_data: IconPath = {
                    "d": path_elem.get("d", ""),
                    "fill": path_elem.get("fill", "none"),  # Lucide typically uses fill="none"
                    "fillRule": path_elem.get("fill-rule"),
                    "clipRule": path_elem.get("clip-rule"),
                }
                # Add stroke attributes from parent SVG if not on path
                if "stroke" not in path_elem.attrib:
                    path_data["stroke"] = stroke  # type: ignore
                    path_data["strokeWidth"] = stroke_width  # type: ignore
                    path_data["strokeLinecap"] = stroke_linecap  # type: ignore
                    path_data["strokeLinejoin"] = stroke_linejoin  # type: ignore
                else:
                    path_data["stroke"] = path_elem.get("stroke")  # type: ignore
                    path_data["strokeWidth"] = path_elem.get("stroke-width")  # type: ignore
                    path_data["strokeLinecap"] = path_elem.get("stroke-linecap")  # type: ignore
                    path_data["strokeLinejoin"] = path_elem.get("stroke-linejoin")  # type: ignore

                # Remove None values
                path_data = {k: v for k, v in path_data.items() if v is not None}  # type: ignore
                paths.append(path_data)  # type: ignore

            return IconData(viewBox=viewBox, paths=paths)

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
