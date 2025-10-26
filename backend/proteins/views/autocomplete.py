from __future__ import annotations

from dal import autocomplete
from django_tomselect.autocompletes import AutocompleteModelView

from ..models import Filter, Lineage, Protein, State

# from django.contrib.postgres.search import TrigramSimilarity


class ProteinAutocomplete(autocomplete.Select2QuerySetView):
    def get_results(self, context):
        """Return data for the 'results' key of the response."""
        return [
            {
                "id": result.slug,
                # 'slug': result.slug,
                "text": result.name,
            }
            for result in context["object_list"]
        ]

    def get_queryset(self):
        if self.request.GET.get("type", "") == "spectra":
            qs = Protein.objects.with_spectra()
        else:
            qs = Protein.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs


class LineageAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        # if not self.request.user.is_authenticated:
        #     return State.objects.none()
        qs = Lineage.objects.all().prefetch_related("protein").order_by("protein__name")
        if self.q:
            qs = qs.filter(protein__name__icontains=self.q)
        return qs


class StateAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        # if not self.request.user.is_authenticated:
        #     return State.objects.none()
        qs = State.objects.all().order_by("protein__name")
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


# Tom-Select autocomplete views
class ProteinTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete view for Protein model."""

    model = Protein
    search_lookups = ["name__icontains", "slug__icontains"]
    ordering = ["name"]
    page_size = 20

    def create_result_dict(self, result):
        return {
            "id": result.slug,
            "text": result.name,
        }

    def get_queryset(self):
        qs = super().get_queryset()
        if self.request.GET.get("type") == "spectra":
            qs = qs.filter(default_state__spectra__isnull=False).distinct()
        return qs


class LineageTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete view for Lineage model."""

    model = Lineage
    search_lookups = ["protein__name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        return {
            "id": result.pk,
            "text": result.protein.name,
        }


class StateTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete view for State model."""

    model = State
    search_lookups = ["protein__name__icontains", "name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        return {
            "id": result.pk,
            "text": str(result),
        }


class FilterTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete view for Filter model."""

    model = Filter
    search_lookups = ["name__icontains", "part__icontains"]
    ordering = ["name"]
    page_size = 20
