from django.contrib.postgres.search import TrigramSimilarity
from django.shortcuts import render, redirect
from django.db.models import Count, Prefetch
from ..filters import ProteinFilter
from ..models import Protein, State
import json


def protein_search(request):
    ''' renders html for protein search page  '''

    if request.GET:
        stateprefetch = Prefetch('states', queryset=State.objects.order_by('-is_dark', 'em_max'))
        f = ProteinFilter(
            request.GET,
            queryset=Protein.visible.annotate(nstates=Count('states')).select_related('default_state')
            .prefetch_related(stateprefetch).order_by('default_state__em_max'))

        # if no hits, but name was provided... try trigram search
        if len(f.qs) == 0:
            name = None
            if 'name__icontains' in f.form.data:
                name = f.form.data['name__icontains']
            elif 'name__iexact' in f.form.data:
                name = f.form.data['name__iexact']
            if name:
                f.recs = Protein.objects.annotate(
                    similarity=TrigramSimilarity('name', name)).filter(
                    similarity__gt=0.2).order_by('-similarity')
        if len(f.qs) == 1:
            return redirect(f.qs.first())
    else:
        f = ProteinFilter(request.GET, queryset=Protein.objects.none())
    return render(request, 'proteins/protein_search.html',
                  {'filter': f,
                   'filter_fields': json.dumps({k: list(set(v)) for k, v in f.Meta.fields.items()}),
                   'filter_operators': json.dumps(f.Meta.operators),
                   'filter_labels': json.dumps(f.Meta.labels),
                   })
