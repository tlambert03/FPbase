from django.core.cache import cache
from django.db.models import Case, F, Value, When
from django.db.models.functions import Concat
from django.http import JsonResponse
from django.shortcuts import render

from fpbase.celery import app
from fpbase.util import is_ajax
from proteins.models import FluorState
from proteins.tasks import calc_fret


def fret_chart(request):
    """renders html for protein spectra page"""
    template = "fret.html"

    if is_ajax(request):
        forster_list = cache.get("forster_list_v2")
        if not forster_list:
            job = cache.get("calc_fret_job_v2")
            if cache.get("calc_fret_job_v2"):
                result = app.AsyncResult(job)
                if result.ready():
                    forster_list = result.get()
                    cache.set("forster_list_v2", forster_list, 60 * 60 * 24)
                    cache.delete("calc_fret_job_v2")
            else:
                job = calc_fret.delay()
                if job:
                    cache.set("calc_fret_job_v2", job.id)
        return JsonResponse({"data": forster_list})

    # Donors need ext_coeff, qy, and an excitation/absorption spectrum.
    # Acceptors only need ext_coeff + an absorption (or excitation) spectrum —
    # dark quenchers like BHQ-10 are legitimate acceptors with no qy.
    base = (
        FluorState.objects.exclude(ext_coeff=None)
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

    def _to_slug_dicts(qs):
        return [
            {
                "slug": x["slug"],
                "category": x["spectra__category"],
                "subtype": x["spectra__subtype"],
                "name": x["display_name"],
            }
            for x in qs
        ]

    donor_slugs = _to_slug_dicts(base.exclude(qy=None))
    acceptor_slugs = _to_slug_dicts(base)

    return render(
        request,
        template,
        {"donor_slugs": donor_slugs, "acceptor_slugs": acceptor_slugs},
    )
