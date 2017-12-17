from .models import Author, Reference
from django.views.generic import DetailView

# Create your views here.


class AuthorDetailView(DetailView):
    ''' renders html for single author page  '''
    model = Author


class ReferenceDetailView(DetailView):
    ''' renders html for single reference page  '''
    model = Reference
