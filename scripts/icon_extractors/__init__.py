"""Base interface for icon extraction from different libraries.

This module defines the standard format for icon data and provides a base
class for icon extractors.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
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


class IconData(TypedDict):
    """Standard icon data format for FPbase.

    This format is library-agnostic and allows easy switching between
    icon libraries.
    """

    viewBox: str  # SVG viewBox attribute (e.g., "0 0 24 24")
    paths: list[IconPath]  # List of SVG paths


class IconExtractor(ABC):
    """Base class for icon extractors.

    Each icon library (FontAwesome, Lucide, etc.) should implement this
    interface to extract icons into the standard IconData format.
    """

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
