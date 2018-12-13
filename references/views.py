from .models import Author, Reference
from django.views.generic import DetailView, ListView
from dal import autocomplete
from django.http import Http404

# Create your views here.


class AuthorDetailView(DetailView):
    ''' renders html for single author page  '''
    queryset = Author.objects.all().prefetch_related('publications', 'publications__authors', 'publications__primary_proteins')


class ReferenceListView(ListView):
    ''' renders html for single reference page  '''
    queryset = Reference.objects.all().prefetch_related('authors', 'proteins', 'primary_proteins')


class ReferenceDetailView(DetailView):
    ''' renders html for single reference page  '''
    queryset = Reference.objects.all().prefetch_related('authors')

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
            except queryset.model.DoesNotExist:
                raise Http404('No reference found matching this query')
            return obj


class ReferenceAutocomplete(autocomplete.Select2QuerySetView):

    def get_results(self, context):
            """Return data for the 'results' key of the response."""
            return [
                {
                    'id': result.doi,
                    # 'slug': result.slug,
                    'text': result.citation,
                } for result in context['object_list']
            ]

    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        if not self.request.user.is_authenticated:
            return Reference.objects.none()
        qs = Reference.objects.all()
        if self.q:
            qs = qs.filter(doi__icontains=self.q)
        return qs
