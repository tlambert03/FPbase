from django import template
from django.template.loader import render_to_string

register = template.Library()


@register.simple_tag(takes_context=True)
def collection_remove_button(context, target):
    user = context['request'].user
    collection = context['proteincollection']

    # do nothing when user isn't authenticated
    if not user.is_authenticated:
        return ''

    if not collection.owner == user:
        return ''

    return render_to_string(
        'proteins/collection_remove_button.html', {
            'target_protein': target.id,
            'target_collection': collection.id,
        }
    )
