from __future__ import annotations

from pathlib import Path

from django import template
from django.utils.html import escape
from django.utils.safestring import mark_safe

register = template.Library()

_ICON_CACHE: dict[str, str] = {}

ICON_DIR = Path(__file__).parent.parent / "static" / "icons"
if not ICON_DIR.exists():
    raise RuntimeError(
        f"Icon directory does not exist: {ICON_DIR}. Please run 'python scripts/extract_fa_icons.py' to generate it."
    )


ICON_KEYS = {svg.stem for svg in ICON_DIR.glob("*.svg")}


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

    if name not in _ICON_CACHE:
        path = ICON_DIR / f"{name}.svg"
        try:
            _ICON_CACHE[name] = path.read_text(encoding="utf8")
        except FileNotFoundError:
            raise ValueError(f"Icon '{name}' not found in FPbase icon vocabulary.") from None
    icon_svg = _ICON_CACHE[name]

    # Build the class attribute (escape to prevent XSS)
    # Build the style attribute (escape to prevent XSS)
    class_ += (class_ + " svg-inline-icon").strip()
    class_attr = f' class="{escape(class_)}"'
    style_attr = f' style="{escape(style)}"' if style else ""

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

    svg = icon_svg.replace("<svg", f"<svg {class_attr}{style_attr}{attrs_str}", 1)
    return mark_safe(svg)
