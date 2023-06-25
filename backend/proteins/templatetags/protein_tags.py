from os.path import splitext

from django import template
from django.template.loader import render_to_string
from django.templatetags.static import static
from django.utils.safestring import mark_safe

register = template.Library()


@register.simple_tag
def webp_picture(name, classed="", alt=""):
    webpname = f"{splitext(name)[0]}.webp"
    return mark_safe(
        f'<picture><source srcset="{static(webpname)}" type="image/webp">'
        f'<img src="{static(name)}" class="{classed}" alt="{alt}"></picture>'
    )


@register.simple_tag(takes_context=True)
def collection_remove_button(context, target):
    user = context["request"].user
    if not user.is_authenticated:
        return ""

    collection = context["proteincollection"]

    return (
        ""
        if collection.owner != user
        else render_to_string(
            "proteins/collection_remove_button.html",
            {"target_protein": target.id, "target_collection": collection.id},
        )
    )


@register.simple_tag(takes_context=True)
def flag_object(context, target):
    user = context["request"].user

    target_model = ".".join((target._meta.app_label, target._meta.object_name))

    return render_to_string(
        "proteins/object_flag_button.html",
        {
            "request": context["request"],
            "target_id": target.id,
            "target_model": target_model,
            "status": target.status,
            "userauth": user.is_authenticated,
        },
    )
