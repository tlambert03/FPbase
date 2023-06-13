import json

from django.contrib.postgres.search import TrigramSimilarity
from django.db.models import Count, Prefetch, Q
from django.shortcuts import redirect, render

from proteins.util.helpers import getprot
from references.models import Author, Reference

from ..filters import ProteinFilter
from ..models import Organism, Protein, State


def protein_search(request):
    """renders html for protein search page"""

    if request.GET:
        if set(request.GET.keys()) == {"q"}:
            query = request.GET.get("q").strip()
            page = None
            try:
                page = getprot(query, visible=True)
                return redirect(page)
            except Protein.DoesNotExist:
                pass
            try:
                page = Protein.visible.get(
                    Q(genbank__iexact=query)
                    | Q(uniprot__iexact=query)
                    | Q(pdb__contains=[query.upper()])
                    | Q(uuid__iexact=query.upper())
                )
                return redirect(page)
            except Protein.DoesNotExist:
                pass
            try:
                page = Author.objects.filter(family__iexact=query).annotate(nr=Count("publications")).order_by("-nr")
                if page:
                    return redirect(page.first())
            except Author.DoesNotExist:
                pass
            try:
                page = Reference.objects.get(doi=query.lower())
                return redirect(page)
            except Reference.DoesNotExist:
                pass
            if len(query) > 5:
                try:
                    page = (
                        Organism.objects.filter(scientific_name__istartswith=query)
                        .annotate(np=Count("proteins"))
                        .order_by("-np")
                    )
                    if page.exists():
                        return redirect(page.first())
                except Organism.DoesNotExist:
                    pass

            request.GET._mutable = True
            request.GET["name__icontains"] = query
            del request.GET["q"]
            return redirect("/search/?name__iexact=" + query)

        stateprefetch = Prefetch("states", queryset=State.objects.order_by("-is_dark", "em_max"))
        f = ProteinFilter(
            request.GET,
            queryset=Protein.visible.annotate(nstates=Count("states"))
            .select_related("default_state")
            .prefetch_related(stateprefetch)
            .order_by("default_state__em_max"),
        )

        # if no hits, but name was provided... try trigram search
        if len(f.qs) == 0:
            name = None
            if "name__icontains" in f.form.data:
                name = f.form.data["name__icontains"]
            elif "name__iexact" in f.form.data:
                name = f.form.data["name__iexact"]
            if name:
                f.recs = (
                    Protein.visible.annotate(similarity=TrigramSimilarity("name", name))
                    .filter(similarity__gt=0.2)
                    .order_by("-similarity")
                )
        if len(f.qs) == 1:
            return redirect(f.qs.first())
    else:
        f = ProteinFilter(request.GET, queryset=Protein.visible.none())
    return render(
        request,
        "proteins/protein_search.html",
        {
            "filter": f,
            "filter_fields": json.dumps(f.Meta.form_fields),
            "filter_operators": json.dumps(f.Meta.operators),
            "filter_labels": json.dumps(f.Meta.labels),
        },
    )
