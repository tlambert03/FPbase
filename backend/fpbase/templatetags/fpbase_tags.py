from __future__ import annotations

from django import template
from django.utils.safestring import mark_safe

register = template.Library()


@register.simple_tag
def icon(name, fa_style="fas", class_="", css_style="", **attrs):
    """Render a FontAwesome icon.

    Parameters
    ----------
    name : str
        The icon name (without 'fa-' prefix). E.g., 'heart', 'external-link-alt'
    fa_style : str, optional
        The FontAwesome style prefix ('fas', 'far', 'fab', 'fa'), by default 'fas'
    class_ : str, optional
        Additional CSS classes to apply, by default ''
    css_style : str, optional
        Inline CSS styles, by default ''
    **attrs
        Additional HTML attributes (e.g., title, data-toggle)

    Returns
    -------
    str
        Safe HTML string for the icon

    Examples
    --------
    >>> {% icon "heart" %}
    <i class="fas fa-heart"></i>

    >>> {% icon "heart" fa_style="far" %}
    <i class="far fa-heart"></i>

    >>> {% icon "heart" class_="mr-2 text-info" %}
    <i class="fas fa-heart mr-2 text-info"></i>

    >>> {% icon "external-link-alt" class_="ml-2" css_style="font-size: 0.8rem;" %}
    <i class="fas fa-external-link-alt ml-2" style="font-size: 0.8rem;"></i>
    """
    # Build the class list
    classes = [fa_style, f"fa-{name}"]
    if class_:
        classes.append(class_)

    # Build the attributes string
    attrs_str = ""
    if css_style:
        attrs_str += f' style="{css_style}"'

    for key, value in attrs.items():
        # Convert underscores to hyphens for data attributes (data_toggle -> data-toggle)
        attr_name = key.replace("_", "-")
        # Handle boolean attributes
        if isinstance(value, bool):
            if value:
                attrs_str += f" {attr_name}"
        else:
            attrs_str += f' {attr_name}="{value}"'

    return mark_safe(f'<i class="{" ".join(classes)}"{attrs_str}></i>')
