from dal import autocomplete
from dal.widgets import QuerySetSelectMixin
from dal_select2.widgets import Select2WidgetMixin
from django import forms


class Select2QuerySetView(autocomplete.Select2QuerySetView):
    pass


class ModelSelect2(QuerySetSelectMixin, Select2WidgetMixin, forms.Select):
    """Select widget for QuerySet choices and Select2."""


class ModelSelect2Multiple(QuerySetSelectMixin, Select2WidgetMixin, forms.SelectMultiple):
    """SelectMultiple widget for QuerySet choices and Select2."""
