from django.core.cache import cache
from django.http import HttpRequest
from django.urls import reverse
from django.utils.cache import get_cache_key


def get_view_cache_key(view_name, args=None, namespace=None, key_prefix=None, request=None):
    """
    This function allows you to invalidate any view-level cache.
    view_name: view function you wish to invalidate or it's named url pattern
    args: any arguments passed to the view function
    namespace: optional, if an application namespace is needed
    key prefix: for the @cache_page decorator for the function (if any)
    """
    # create a proxy request object to access cache
    if request:
        req = request
    else:
        req = HttpRequest()
        req.META = {"HTTP_HOST": "127.0.0.1:8000", "SERVER_PORT": 8000}

    if args is None:
        args = []

    # Loookup the req path:
    if namespace:
        view_name = namespace + ":" + view_name
    req.path = reverse(view_name, args=args)

    # get cache key, expire if the cached item exists:
    return get_cache_key(req, key_prefix=key_prefix)


def clear_view_cache(*args, **kwargs):
    cache_key = get_view_cache_key(*args, **kwargs)
    key_deleted = False
    if cache_key:
        key_deleted = cache.delete(cache_key)
    return key_deleted


def uncache_protein_page(slug, request):
    clear_view_cache("proteins:protein-detail", args=[slug], request=request)


def show_queries():
    import logging

    logger = logging.getLogger("django.db.backends")
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())


def is_ajax(request):
    # https://stackoverflow.com/a/70419609
    return request.headers.get("x-requested-with") == "XMLHttpRequest"
