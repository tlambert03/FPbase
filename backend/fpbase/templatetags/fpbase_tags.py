from __future__ import annotations

import json
from pathlib import Path

from django import template
from django.conf import settings
from django.utils.html import escape
from django.utils.safestring import mark_safe

register = template.Library()

# Load icon data from JSON file based on configured library
LIBRARY = "fontawesome"
icon_file = Path(settings.APPS_DIR) / "static" / "icons" / f"{LIBRARY}.json"
if not icon_file.exists():
    raise FileNotFoundError(
        f"Icon data file not found: {icon_file}. Run 'python scripts/extract_icons.py {LIBRARY}' to generate it."
    )

ICON_LIBRARY = json.loads(icon_file.read_text())
ICON_DATA = ICON_LIBRARY.get("icons", {})
LIBRARY_DEFAULTS = ICON_LIBRARY.get("defaults", {})
ICON_SCALE = ICON_LIBRARY.get("scale", "1em")


# Available icon names (for validation and error messages)
AVAILABLE_ICONS = {
    "info",
    "info-circle",
    "warning",
    "alert",
    "help",
    "question",
    "close",
    "remove",
    "menu",
    "grid",
    "search",
    "filter",
    "view",
    "settings",
    "edit",
    "delete",
    "trash",
    "undo",
    "check",
    "success",
    "selected",
    "unselected",
    "heart",
    "heart-outline",
    "add",
    "add-item",
    "download",
    "upload",
    "share",
    "share-square",
    "link",
    "external-link",
    "exchange",
    "book",
    "collection",
    "quote",
    "photo",
    "chart",
    "table",
    "flag",
    "flag-outline",
    "clock",
    "spinner",
    "lightbulb",
    "sun",
    "email",
    "wrench",
    "keyboard",
    "google",
    "twitter",
    "orcid",
}


@register.simple_tag
def icon(name, class_="", style="", **attrs):
    """Render an icon using FPbase's icon vocabulary.

    This renders inline SVG icons from the configured icon library, eliminating
    the need for external icon font CDNs.

    Parameters
    ----------
    name : str
        The FPbase icon name (e.g., 'info', 'warning', 'external-link')
    class_ : str, optional
        Additional CSS classes to apply to the SVG, by default ''
    style : str, optional
        Inline CSS styles, by default ''
    **attrs
        Additional HTML attributes (e.g., title, data_toggle, aria_hidden)

    Raises
    ------
    ValueError
        If the icon name is not found in the icon data

    Returns
    -------
    str
        Safe HTML string containing the SVG icon
    """
    if name not in AVAILABLE_ICONS:
        raise ValueError(
            f"Icon '{name}' not found in FPbase icon vocabulary. Available icons: {', '.join(sorted(AVAILABLE_ICONS))}"
        )

    # Load icon data
    if name not in ICON_DATA:
        raise ValueError(
            f"Icon '{name}' not available in {LIBRARY} library. "
            "Try switching to a different library or use an alternative icon."
        )

    icon_svg = ICON_DATA[name]
    viewBox = icon_svg["viewBox"]
    paths = icon_svg["paths"]

    # Build the class attribute (escape to prevent XSS)
    class_attr = ""
    if class_:
        class_attr = f' class="{escape(class_)}"'

    # Build the style attribute (escape to prevent XSS)
    style_attr = ""
    if style:
        style_attr = f' style="{escape(style)}"'

    # Build additional attributes (with XSS protection)
    attrs_str = ""
    for key, value in attrs.items():
        # Convert underscores to hyphens for data attributes (data_toggle -> data-toggle)
        attr_name = key.replace("_", "-")
        # Handle boolean attributes
        if isinstance(value, bool):
            if value:
                attrs_str += f" {attr_name}"
        else:
            attrs_str += f' {attr_name}="{escape(str(value))}"'

    # Merge library defaults with icon-specific attributes
    # Icon-specific attributes override defaults
    merged_attrs = {
        **LIBRARY_DEFAULTS,
        **{
            k: v
            for k, v in icon_svg.items()
            if k in ("stroke", "strokeWidth", "strokeLinecap", "strokeLinejoin", "fill")
        },
    }

    # Extract SVG-level attributes (for stroke-based icons like Lucide)
    svg_attrs = []
    for attr_name in ("stroke", "strokeWidth", "strokeLinecap", "strokeLinejoin", "fill"):
        if attr_name in merged_attrs:
            # Convert camelCase to kebab-case
            kebab_name = (
                "stroke-width"
                if attr_name == "strokeWidth"
                else "stroke-linecap"
                if attr_name == "strokeLinecap"
                else "stroke-linejoin"
                if attr_name == "strokeLinejoin"
                else attr_name
            )
            svg_attrs.append(f'{kebab_name}="{escape(str(merged_attrs[attr_name]))}"')

    svg_attrs_str = " " + " ".join(svg_attrs) if svg_attrs else ""

    # Build path elements
    path_elements = []
    for path in paths:
        path_attrs = []
        for key, value in path.items():
            # Convert camelCase to kebab-case for SVG attributes
            if key == "fillRule":
                attr_name = "fill-rule"
            elif key == "clipRule":
                attr_name = "clip-rule"
            elif key == "strokeWidth":
                attr_name = "stroke-width"
            elif key == "strokeLinecap":
                attr_name = "stroke-linecap"
            elif key == "strokeLinejoin":
                attr_name = "stroke-linejoin"
            else:
                attr_name = key

            path_attrs.append(f'{attr_name}="{escape(str(value))}"')

        path_elements.append(f"<path {' '.join(path_attrs)}/>")

    svg = (
        f"<svg{class_attr}{style_attr} "
        f'xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="{viewBox}" '
        f'width="{ICON_SCALE}" '
        f'height="{ICON_SCALE}"{svg_attrs_str} '
        f'aria-hidden="true"{attrs_str}>'
        f"{''.join(path_elements)}"
        f"</svg>"
    )

    return mark_safe(svg)
