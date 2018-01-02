from .models import Author, Reference
from django.views.generic import DetailView

# Create your views here.


class AuthorDetailView(DetailView):
    ''' renders html for single author page  '''
    queryset = Author.objects.all().prefetch_related('publications', 'publications__authors', 'publications__primary_proteins')


class ReferenceDetailView(DetailView):
    ''' renders html for single reference page  '''
    queryset = Reference.objects.all().prefetch_related('authors')
