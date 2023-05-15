from django.core.cache import cache
from django.http import JsonResponse
from django.shortcuts import render

from fpbase.celery import app
from fpbase.util import is_ajax

from ..models import Dye, State
from ..tasks import calc_fret


def fret_chart(request):
    """renders html for protein spectra page"""
    template = "fret.html"

    if is_ajax(request):
        forster_list = cache.get("forster_list")
        if not forster_list:
            job = cache.get("calc_fret_job")
            if cache.get("calc_fret_job"):
                result = app.AsyncResult(job)
                if result.ready():
                    forster_list = result.get()
                    cache.set("forster_list", forster_list, 60 * 60 * 24)
                    cache.delete("calc_fret_job")
            else:
                job = calc_fret.delay()
                if job:
                    cache.set("calc_fret_job", job.id)
        return JsonResponse({"data": forster_list})

    slugs = (
        State.objects.exclude(ext_coeff=None)
        .exclude(qy=None)
        .filter(spectra__subtype__in=("ex", "ab"))
        .values("slug", "name", "protein__name", "spectra__category", "spectra__subtype")
    )

    slugs = [
        {
            "slug": x["slug"],
            "category": x["spectra__category"],
            "subtype": x["spectra__subtype"],
            "name": x["protein__name"] + (f" ({x['name']})" if x["name"] != "default" else ""),
        }
        for x in slugs
    ]

    good_dyes = (
        Dye.objects.exclude(ext_coeff=None).exclude(qy=None).filter(spectra__subtype__in=("ex", "ab"))
    ).values("slug", "name", "spectra__category", "spectra__subtype")

    slugs += [
        {
            "slug": x["slug"],
            "category": x["spectra__category"],
            "subtype": x["spectra__subtype"],
            "name": x["name"],
        }
        for x in good_dyes
    ]

    slugs.sort(key=lambda x: x["name"])
    return render(request, template, {"probeslugs": slugs})
