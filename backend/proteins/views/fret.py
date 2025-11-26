from django.core.cache import cache
from django.db.models import Case, F, Value, When
from django.db.models.functions import Concat
from django.http import JsonResponse
from django.shortcuts import render

from fpbase.celery import app
from fpbase.util import is_ajax

from ..models import Fluorophore
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

    # Query all fluorophores (States + DyeStates) with required properties
    # Build display name and sort in the database
    fluorophores = (
        Fluorophore.objects.exclude(ext_coeff=None)
        .exclude(qy=None)
        .filter(spectra__subtype__in=("ex", "ab"))
        .annotate(
            display_name=Case(
                When(name="default", then=F("owner_name")),
                default=Concat(F("owner_name"), Value(" ("), F("name"), Value(")")),
            )
        )
        .order_by("display_name")
        .values("slug", "display_name", "spectra__category", "spectra__subtype")
    )

    slugs = [
        {
            "slug": x["slug"],
            "category": x["spectra__category"],
            "subtype": x["spectra__subtype"],
            "name": x["display_name"],
        }
        for x in fluorophores
    ]

    return render(request, template, {"probeslugs": slugs})
