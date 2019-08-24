from django.core.cache import cache
from django.http import JsonResponse
from django.shortcuts import render

from fpbase.celery import app

from ..models import State
from ..tasks import calc_fret


def fret_chart(request):
    """ renders html for protein spectra page  """
    template = "fret.html"

    if request.is_ajax():
        L = cache.get("forster_list")
        if not L:
            job = cache.get("calc_fret_job")
            if cache.get("calc_fret_job"):
                result = app.AsyncResult(job)
                if result.ready():
                    L = result.get()
                    cache.set("forster_list", L, 60 * 60 * 24)
                    cache.delete("calc_fret_job")
            else:
                job = calc_fret.delay()
                if job:
                    cache.set("calc_fret_job", job.id)
        return JsonResponse({"data": L})

    slugs = State.objects.filter(spectra__subtype__in=("ex", "ab")).values(
        "slug", "name", "protein__name", "spectra__category", "spectra__subtype"
    )

    slugs = [
        {
            "slug": x["slug"],
            "category": x["spectra__category"],
            "subtype": x["spectra__subtype"],
            "name": x["protein__name"]
            + (f" ({x['name']})" if x["name"] != "default" else ""),
        }
        for x in slugs
    ]
    return render(request, template, {"probeslugs": slugs})
