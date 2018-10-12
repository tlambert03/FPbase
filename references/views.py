from .models import Author, Reference
from django.views.generic import DetailView, ListView
from dal import autocomplete

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
