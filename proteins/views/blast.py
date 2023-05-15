from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import ensure_csrf_cookie

from fpbase.util import is_ajax

from ..util.blast import blast


@ensure_csrf_cookie
def blast_view(request):
    if not is_ajax(request):
        return render(request, "proteins/blast.html")

    seq = request.POST.get("query")
    binary = request.POST.get("binary", "blastp")
    assert binary in ("blastx", "blastp")
    if seq:
        try:
            result = blast(seq, binary)
            return JsonResponse({"status": 200, "blastResult": result})
        except Exception as e:
            return JsonResponse({"status": 500, "error": f"BLAST error: {e}"})
    return JsonResponse({"status": 204, "blastResult": []})
