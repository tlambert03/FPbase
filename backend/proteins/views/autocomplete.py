from tomcomplete import Select2QuerySetView

from ..models import Filter, Lineage, Protein, State

# from django.contrib.postgres.search import TrigramSimilarity


class ProteinAutocomplete(Select2QuerySetView):
    def get_results(self, context):
        """Return data for the 'results' key of the response."""
        return [{"id": result.slug, "text": result.name} for result in context["object_list"]]

    def get_queryset(self):
        if self.request.GET.get("type", "") == "spectra":
            qs = Protein.objects.with_spectra()
        else:
            qs = Protein.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs


class LineageAutocomplete(Select2QuerySetView):
    def get_queryset(self):
        qs = Lineage.objects.all().prefetch_related("protein").order_by("protein__name")
        if self.q:
            qs = qs.filter(protein__name__icontains=self.q)
        return qs


class StateAutocomplete(Select2QuerySetView):
    def get_queryset(self):
        qs = State.objects.all().order_by("protein__name")
        if self.q:
            qs = qs.filter(protein__name__icontains=self.q)
        return qs


class FilterAutocomplete(Select2QuerySetView):
    def get_queryset(self):
        qs = Filter.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs
