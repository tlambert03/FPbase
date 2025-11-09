from __future__ import annotations

from django import template
from django.utils.html import escape
from django.utils.safestring import mark_safe

register = template.Library()

# FPbase icon vocabulary mapped to FontAwesome
# This allows us to swap icon libraries in the future without changing templates
ICON_MAP = {
    # UI & Navigation
    "info": ("fas", "info-circle"),
    "info-i": ("fas", "info"),
    "warning": ("fas", "exclamation-circle"),
    "alert": ("fas", "exclamation-triangle"),
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


@register.simple_tag
def icon(name, class_="", style="", **attrs):
    """Render an icon using FPbase's icon vocabulary.

    This abstraction allows FPbase to use semantic icon names that can be mapped
    to any icon library (currently FontAwesome, but could be Lucide, etc. in the future).

    Parameters
    ----------
    name : str
        The FPbase icon name (e.g., 'info', 'warning', 'external-link')
        See ICON_MAP for available icons.
    class_ : str, optional
        Additional CSS classes to apply, by default ''
    style : str, optional
        Inline CSS styles, by default ''
    **attrs
        Additional HTML attributes (e.g., title, data_toggle, aria_hidden)

    Raises
    ------
    ValueError
        If the icon name is not found in ICON_MAP
    """
    if name not in ICON_MAP:
        raise ValueError(
            f"Icon '{name}' not found in FPbase icon vocabulary. Available icons: {', '.join(sorted(ICON_MAP.keys()))}"
        )

    # Get the icon library implementation (currently FontAwesome)
    fa_style, fa_icon_name = ICON_MAP[name]

    # Build the class list (escape class_ to prevent XSS)
    classes = [fa_style, f"fa-{fa_icon_name}"]
    if class_:
        classes.append(escape(class_))

    # Build the attributes string (with XSS protection via escaping)
    attrs_str = ""
    if style:
        attrs_str += f' style="{escape(style)}"'

    for key, value in attrs.items():
        # Convert underscores to hyphens for data attributes (data_toggle -> data-toggle)
        attr_name = key.replace("_", "-")
        # Handle boolean attributes
        if isinstance(value, bool):
            if value:
                attrs_str += f" {attr_name}"
        else:
            attrs_str += f' {attr_name}="{escape(str(value))}"'

    return mark_safe(f'<i class="{" ".join(classes)}"{attrs_str}></i>')
