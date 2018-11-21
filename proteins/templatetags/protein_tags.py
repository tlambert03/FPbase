from django import template
from django.conf import settings
from django.template.loader import render_to_string
from django.contrib.staticfiles.storage import staticfiles_storage
from webpack_loader.templatetags.webpack_loader import render_bundle
from webpack_loader.utils import _get_bundle

register = template.Library()


@register.simple_tag
def custom_static(bundle_name, extension=None, config='DEFAULT', attrs=''):
    if settings.DEBUG:
        return render_bundle(bundle_name, extension=extension, config=config, attrs=attrs)

    bundle = _get_bundle(bundle_name, extension, config)
    tags = []
    for chunk in bundle:
        if chunk['name'].endswith(('.js', '.js.gz')):
            tags.append((
                '<script type="text/javascript" src="{0}" {1}></script>'
            ).format(staticfiles_storage.url(chunk['name']), attrs))
        elif chunk['name'].endswith(('.css', '.css.gz')):
            tags.append((
                '<link type="text/css" href="{0}" rel="stylesheet" {1}/>'
            ).format(staticfiles_storage.url(chunk['name']), attrs))
    return "\n".join(tags)


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


@register.simple_tag(takes_context=True)
def flag_object(context, target):
    user = context['request'].user

    target_model = '.'.join((target._meta.app_label, target._meta.object_name))

    return render_to_string(
        'proteins/object_flag_button.html', {
            'request': context['request'],
            'target_id': target.id,
            'target_model': target_model,
            'status': target.status,
            'userauth': user.is_authenticated
        }
    )
