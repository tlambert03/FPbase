from __future__ import annotations

from django_tomselect.autocompletes import AutocompleteModelView

from ..models import Filter, Lineage, Protein, State


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
