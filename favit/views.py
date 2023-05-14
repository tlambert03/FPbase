from django.apps import apps
from django.contrib.auth.decorators import login_required
from django.http import HttpResponseBadRequest, HttpResponseNotAllowed, JsonResponse

from fpbase.util import is_ajax, uncache_protein_page

from .models import Favorite


@login_required
def add_or_remove(request):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])
    user = request.user
    try:
        app_model = request.POST["target_model"]
        obj_id = int(request.POST["target_object_id"])
    except (KeyError, ValueError):
        return HttpResponseBadRequest()
    fav = Favorite.objects.get_favorite(user, obj_id, model=app_model)
    if fav is None:
        Favorite.objects.create(user, obj_id, app_model)
        status = "added"
    else:
        fav.delete()
        status = "deleted"
    response = {
        "status": status,
        "fav_count": Favorite.objects.for_object(obj_id, app_model).count(),
    }

    if str(app_model) == "proteins.Protein":
        try:
            slug = apps.get_model(app_model).objects.get(id=obj_id).slug
            uncache_protein_page(slug, request)
        except Exception as e:
            import logging

            logger = logging.getLogger(__name__)
            logger.error(f"failed to uncache protein: {e}")

    return JsonResponse(response)


@login_required
def remove(request):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])

    user = request.user

    try:
        app_model = request.POST["target_model"]
        obj_id = int(request.POST["target_object_id"])
    except (KeyError, ValueError):
        return HttpResponseBadRequest()

    Favorite.objects.get_favorite(user, obj_id, model=app_model).delete()
    status = "deleted"

    response = {"status": status}

    return JsonResponse(response)
