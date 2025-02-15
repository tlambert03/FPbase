import reversion
from dal import autocomplete
from django.core.mail import mail_managers
from django.http import Http404, HttpResponseNotAllowed, JsonResponse
from django.utils.html import strip_tags
from django.views.generic import DetailView, ListView

from fpbase.util import is_ajax
from proteins.models import Excerpt
from proteins.util.helpers import link_excerpts

from .models import Author, Reference


class AuthorDetailView(DetailView):
    """renders html for single author page"""

    queryset = Author.objects.all().prefetch_related(
        "publications", "publications__authors", "publications__primary_proteins"
    )


class ReferenceListView(ListView):
    """renders html for single reference page"""

    queryset = Reference.objects.all().prefetch_related("authors", "proteins", "primary_proteins")


class ReferenceDetailView(DetailView):
    """renders html for single reference page"""

    queryset = Reference.objects.all().prefetch_related("authors")

    def get_object(self, queryset=None):
        try:
            return super().get_object(queryset=queryset)
        except (ValueError, Http404):
            # allow for doi to be used in url as well
            if queryset is None:
                queryset = self.get_queryset()
            try:
                doi = self.kwargs.get(self.pk_url_kwarg)
                queryset = queryset.filter(doi=doi.lower())
                obj = queryset.get()
            except queryset.model.DoesNotExist as e:
                raise Http404("No reference found matching this query") from e
            return obj

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        data["excerpts"] = link_excerpts(self.object.excerpts.all())
        return data


class ReferenceAutocomplete(autocomplete.Select2QuerySetView):
    def get_results(self, context):
        """Return data for the 'results' key of the response."""
        return [
            {
                "id": result.doi,
                # 'slug': result.slug,
                "text": result.citation,
            }
            for result in context["object_list"]
        ]

    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        if not self.request.user.is_authenticated:
            return Reference.objects.none()
        qs = Reference.objects.all()
        if self.q:
            qs = qs.filter(doi__icontains=self.q)
        return qs


def add_excerpt(request, pk=None):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])
    try:
        with reversion.create_revision():
            ref = Reference.objects.get(pk=pk)
            content = request.POST.get("excerpt_content")
            if content:
                # P.references.add(ref)
                Excerpt.objects.create(reference=ref, content=strip_tags(content), created_by=request.user)
                if not request.user.is_staff:
                    msg = (
                        f"User: {request.user.username}\nReference: {ref}, {ref.title}"
                        f"\nExcerpt: {strip_tags(content)}\n{request.build_absolute_uri(ref.get_absolute_url())}"
                    )
                    mail_managers("Excerpt Added", msg, fail_silently=True)
                reversion.set_user(request.user)
                reversion.set_comment(f"Excerpt from {ref} added")
        return JsonResponse({"status": "success"})
    except Exception as e:
        return JsonResponse({"status": "failed", "msg": e})
