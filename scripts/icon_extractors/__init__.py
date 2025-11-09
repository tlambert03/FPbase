"""Base interface for icon extraction from different libraries.

This module defines the standard format for icon data and provides a base
class for icon extractors.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Required, TypedDict


class IconPath(TypedDict, total=False):
    """A single SVG path within an icon.

    Supports both fill-based icons (FontAwesome) and stroke-based icons (Lucide).
    All fields except 'd' are optional.
    """

    d: Required[str]  # SVG path data (required)
    fill: str  # Fill color (usually "currentColor" or "none")
    fillRule: str  # Fill rule (e.g., "evenodd")
    clipRule: str  # Clip rule
    stroke: str  # Stroke color (for stroke-based icons)
    strokeWidth: str  # Stroke width
    strokeLinecap: str  # Stroke linecap (e.g., "round")
    strokeLinejoin: str  # Stroke linejoin (e.g., "round")


class IconData(TypedDict, total=False):
    """Standard icon data format for FPbase.

    This format is library-agnostic and allows easy switching between
    icon libraries.
    """

    viewBox: Required[str]  # SVG viewBox attribute (e.g., "0 0 24 24")
    paths: Required[list[IconPath]]  # List of SVG paths

    # Optional SVG-level attributes (used by stroke-based icons like Lucide)
    # These apply to all paths unless overridden at the path level
    stroke: str  # Default stroke color
    strokeWidth: str  # Default stroke width
    strokeLinecap: str  # Default stroke linecap
    strokeLinejoin: str  # Default stroke linejoin
    fill: str  # Default fill color


class IconLibrary(TypedDict, total=False):
    """Top-level structure for icon library JSON files.

    This allows each library to specify its own defaults while keeping
    the template tag implementation library-agnostic.
    """

    icons: Required[dict[str, IconData]]  # Icon name -> icon data mapping
    defaults: dict[str, str]  # Default SVG attributes for this library
    scale: str  # Optional scale override (e.g., "1.1em" for Lucide)


class IconExtractor(ABC):
    """Base class for icon extractors.

    Each icon library (FontAwesome, Lucide, etc.) should implement this
    interface to extract icons into the standard IconData format.
    """

    @classmethod
    def extract_to(cls, output_file: Path) -> None:
        self = cls()
        library_data = self.extract_library()
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(json.dumps(library_data, indent=2))
        print(f"\n✓ Extracted {len(library_data['icons'])} icons to {output_file}")

    @abstractmethod
    def extract_icon(self, icon_name: str) -> IconData | None:
        """Extract a single icon by name.

        Parameters
        ----------
        icon_name : str
            The library-specific icon name (e.g., "external-link-alt" for FA)

        Returns
        -------
        IconData | None
            The extracted icon data in standard format, or None if not found
        """
        pass

    @abstractmethod
    def get_available_icons(self) -> list[str]:
        """Get a list of all available icon names in the library.

        Returns
        -------
        list[str]
            List of icon names available in this library
        """
        pass

    @abstractmethod
    def extract_library(self) -> IconLibrary:
        """Extract all icons for this library.

        This method should extract all icons defined in the library's
        ICON_MAP and return a complete IconLibrary structure ready to
        be written to JSON.

        Returns
        -------
        IconLibrary
            Complete library data including icons, defaults, and optional scale
        """
        pass

    def extract_icons(self, icon_names: list[str]) -> dict[str, IconData]:
        """Extract multiple icons by name.

        Parameters
        ----------
        icon_names : list[str]
            List of icon names to extract

        Returns
        -------
        dict[str, IconData]
            Dictionary mapping icon names to their data
        """
        icons = {}
        for name in icon_names:
            icon_data = self.extract_icon(name)
            if icon_data:
                icons[name] = icon_data
            else:
                print(f"Warning: Icon '{name}' not found in library")
        return icons
