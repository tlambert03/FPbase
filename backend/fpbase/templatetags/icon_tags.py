"""Template tags for rendering SVG icons."""

from __future__ import annotations

from django import template
from django.utils.safestring import mark_safe

from fpbase.icons import get_icon

register = template.Library()


@register.simple_tag
def icon(name: str, class_name: str = "", **attrs) -> str:
    """
    Render an SVG icon.

    Usage:
        {% load icon_tags %}
        {% icon 'trash' %}
        {% icon 'trash' class='text-danger' %}
        {% icon 'edit' class='fa-lg' %}
        {% icon 'spinner' class='fa-spin' aria_label='Loading' %}

    Args:
        name: Icon name (kebab-case like 'trash' or camelCase like 'fasTrash')
        class_name: Optional CSS classes to add to the SVG
        **attrs: Additional HTML attributes (use underscores instead of hyphens)

    Returns:
        SVG markup string
    """
    icon_data = get_icon(name)

    if not icon_data:
        # Return a placeholder or warning in dev mode
        return mark_safe(f'<!-- Icon "{name}" not found -->')

    # Build attributes dict
    svg_attrs = {
        "viewBox": f"0 0 {icon_data['width']} {icon_data['height']}",
        "fill": "currentColor",
        "aria-hidden": "true",
    }

    # Add custom attributes (replace underscores with hyphens for HTML)
    for key, value in attrs.items():
        html_key = key.replace("_", "-")
        svg_attrs[html_key] = value

    # Build attribute string
    attr_string = " ".join(f'{k}="{v}"' for k, v in svg_attrs.items())

    # Add class if provided
    class_attr = f' class="{class_name}"' if class_name else ""

    # Render SVG
    svg = f'<svg {attr_string}{class_attr}><path d="{icon_data["path"]}"/></svg>'

    return mark_safe(svg)


@register.inclusion_tag("fpbase/icon.html")
def icon_button(icon_name: str, text: str = "", btn_class: str = "btn-primary", **attrs) -> dict:
    """
    Render a button with an icon.

    Usage:
        {% icon_button 'trash' 'Delete' btn_class='btn-danger' %}
        {% icon_button 'plus' 'Add New' %}

    Args:
        icon_name: Icon name
        text: Button text (optional)
        btn_class: Bootstrap button class
        **attrs: Additional button attributes

    Returns:
        Context dict for template
    """
    return {
        "icon": get_icon(icon_name),
        "text": text,
        "btn_class": btn_class,
        "attrs": attrs,
    }
