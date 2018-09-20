from django.shortcuts import render
from django.core.cache import cache
from django.http import JsonResponse
from ..models import State
from ..util.helpers import forster_list


def fret_chart(request):
    ''' renders html for protein spectra page  '''
    template = 'fret.html'

    if request.method == 'GET':

        if request.is_ajax():
            L = cache.get('forster_list')
            if not L:
                L = forster_list()
                cache.set('forster_list', L, 60 * 60 * 24 * 7)
            return JsonResponse({'data': L})

        else:
            slugs = State.objects.filter(spectra__subtype='ex').values('slug', 'protein__name', 'spectra__category', 'spectra__subtype')
            slugs = [{'slug': x['slug'], 'category': x['spectra__category'], 'subtype': x['spectra__subtype'], 'name': x['protein__name']} for x in slugs]
            return render(request, template, {
                'probeslugs': slugs}
            )
