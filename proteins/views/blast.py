from django.shortcuts import render
from django.http import JsonResponse
from ..util.blast import blast
from django.views.decorators.csrf import ensure_csrf_cookie


@ensure_csrf_cookie
def blast_view(request):
    if request.is_ajax():
        seq = request.POST.get('query')
        binary = request.POST.get('binary', 'blastp')
        assert binary in ('blastx', 'blastp')
        if seq:
            try:
                return JsonResponse({
                    'status': 200,
                    'blastResult': blast(seq, binary)
                })
            except Exception:
                return JsonResponse({'status': 500})
        return JsonResponse({
            'status': 204,
            'blastResult': []
        })
    return render(request, 'proteins/blast.html')
