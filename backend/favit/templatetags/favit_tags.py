from django import template
from django.contrib.contenttypes.models import ContentType
from django.template.loader import render_to_string

from ..models import Favorite

register = template.Library()


@register.simple_tag(takes_context=True)
def favorite_button(context, target):
    user = context["request"].user

    # moved to the template
    if not user.is_authenticated:
        return render_to_string(
            "favit/button.html",
            {
                "authenticated": False,
                "path": context["request"].path,
                "fav_count": Favorite.objects.for_object(target).count(),
            },
        )

    target_model = ".".join((target._meta.app_label, target._meta.object_name))

    undo = False
    # prepare button to unfave if the user
    # already faved this object
    if Favorite.objects.get_favorite(user, target):
        undo = True

    return render_to_string(
        "favit/button.html",
        {
            "target_model": target_model,
            "target_object_id": target.id,
            "undo": undo,
            "authenticated": True,
            "fav_count": Favorite.objects.for_object(target).count(),
        },
    )


@register.simple_tag(takes_context=True)
def unfave_button(context, target):
    user = context["request"].user

    # do nothing when user isn't authenticated
    if not user.is_authenticated:
        return ""

    if Favorite.objects.get_favorite(user, target) is None:
        return ""

    target_model = ".".join((target._meta.app_label, target._meta.object_name))

    return render_to_string(
        "favit/unfave-button.html",
        {"target_model": target_model, "target_object_id": target.id},
    )


@register.filter
def get_favorite_for(obj, user):
    """
    Get Favorite instance for an object (obj) and a user (user)

    Usage:
    {% with obj|get_favorite_for:user as fav_object %}
        ...
    {% endwith %}
    """

    return Favorite.objects.get_favorite(user, obj)


@register.filter
def favorites_count(obj):
    """
    Usage:

    Given an object `obj` you may show it fav count like this:

    <p>Favorite Count {{ obj|favorites_count }}</p>
    """

    return Favorite.objects.for_object(obj).count()


@register.filter
def get_object(fav):
    ct_id = fav.target_content_type_id
    model = ContentType.objects.get_for_id(ct_id).model_class()
    try:
        return model.objects.get(id=fav.target_object_id)
    except model.DoesNotExist:
        return None


@register.simple_tag
def user_favorites(user, app_model="proteins.Protein"):
    """
    Usage:

    Get all user favorited objects:

        {% with user_favorites <user> as favorite_list %}
            {% for fav_obj in favorite_list %}
                {# do something with fav_obj #}
            {% endfor %}
        {% endwith %}

    or, just favorites from one model:

        {% with user_favorites <user> "app_label.model" as favorite_list %}
            {% for fav_obj in favorite_list %}
                {# do something with fav_obj #}
            {%
        {% endwith %}
    """

    return Favorite.objects.for_user(user, app_model)


@register.simple_tag
def model_favorites(app_model):
    """
    Gets all favorited objects that are instances of a model
    given in module notation.

    Usage:

        {% with model_favorites "app_label.model" as favorite_list %}
            {% for fav_obj in favorite_list %}
                {# do something with fav_obj #}
            {% endfor %}
        {% endwith %}
    """

    return Favorite.objects.for_model(app_model)
