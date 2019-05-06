from django.shortcuts import render
from django.http import JsonResponse
from ..util.blast import blast


def blast_view(request):
    if request.is_ajax():
        seq = request.POST.get('query')
        binary = request.POST.get('binary', 'blastp')
        assert binary in ('blastx', 'blastp')
        if seq:
            return JsonResponse({
                'status': 200,
                'blastResult': blast(seq, binary)
            })
        return JsonResponse({
            'status': 204,
            'blastResult': []
        })
    return render(request, 'proteins/blast.html')
