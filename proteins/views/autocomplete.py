from dal import autocomplete
from ..models import Protein, State, Filter
# from django.contrib.postgres.search import TrigramSimilarity


class ProteinAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        if not self.request.user.is_authenticated:
            return Protein.objects.none()
        qs = Protein.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs


class StateAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        if not self.request.user.is_authenticated:
            return State.objects.none()
        qs = State.objects.all()
        if self.q:
            qs = qs.filter(protein__name__icontains=self.q)
        return qs


class FilterAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        # if not self.request.user.is_authenticated:
        #     return Filter.objects.none()
        qs = Filter.objects.all()
        if self.q:
            # qs = Filter.objects.annotate(
            #     similarity=TrigramSimilarity('part', self.q)) \
            #     .filter(similarity__gt=0.3) \
            #     .order_by('-similarity')
            qs = qs.filter(name__icontains=self.q)
        return qs
