import copy
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, TypedDict, Unpack

from dal import autocomplete
from dal.widgets import WidgetMixin
from dal_select2.widgets import Select2WidgetMixin
from django import forms
from django.urls import reverse


class QuerySetSelectMixin(WidgetMixin):
    """QuerySet support for choices."""

    def filter_choices_to_render(self, selected_choices: Iterable) -> None:
        """Filter out un-selected choices if choices is a QuerySet."""
        try:
            self.choices.queryset = self.choices.queryset.filter(pk__in=[c for c in selected_choices if c])
        except (ValueError, AttributeError):
            # if selected_choices are invalid, do nothing
            pass


class Select2QuerySetView(autocomplete.Select2QuerySetView):
    pass


class AutocompleteKwargs(TypedDict, total=False):
    value_field: str
    label_field: str
    min_length: int
    delay: int
    placeholder: str
    depends_on: str
    create_url: str


class AutocompleteSelect(forms.Select):
    """A Select widget that adds data attributes for JS enhancement.


    Parameters
    ----------
    autocomplete_url : str
        The URL to use to fetch autocomplete results.
    attrs : Mapping | None
        HTML attributes for the widget.
    choices : Sequence[tuple[Any, Any]]
        The choices for the select widget.
    **kwargs : Any
        Autocomplete configuration options:
        - value_field: JSON field for option value (default: 'id')
        - label_field: JSON field for option label (default: 'text')
        - min_length: Minimum characters before search (default: 2)
        - delay: Debounce delay in ms (default: 300)
        - placeholder: Placeholder text for search input
        - depends_on: Field name this depends on
        - create_url: URL for creating new options
    """

    def __init__(
        self,
        url: str | None = None,
        attrs: Mapping | None = None,
        choices: Sequence[tuple[Any, Any]] = (),
        **autocomplete_config: Unpack[AutocompleteKwargs],
    ):
        """Instanciate a widget with a URL and a list of fields to forward."""
        self.autocomplete_url = url
        self.url = url  # Set for DAL's optgroups filtering logic
        self.autocomplete_config = autocomplete_config
        valid_kwargs = set(AutocompleteKwargs.__annotations__)
        invalid = set(autocomplete_config) - valid_kwargs
        if invalid:
            raise TypeError(f"Invalid autocomplete configuration keys: {', '.join(invalid)}")
        super().__init__(attrs, choices)

    def build_attrs(self, base_attrs, extra_attrs=None):
        """Add data-autocomplete-* attributes to the widget."""
        attrs = super().build_attrs(base_attrs, extra_attrs)

        if url := self.autocomplete_url:
            normed = url if "/" in url else reverse(url)
            attrs["data-autocomplete-url"] = normed

        # Add optional configuration as data-autocomplete- attributes
        for config_key in AutocompleteKwargs.__annotations__:
            attr_name = f"data-autocomplete-{config_key.replace('_', '-')}"
            if config_key in self.autocomplete_config:
                attrs[attr_name] = self.autocomplete_config[config_key]

        return attrs

    def optgroups(self, name: str, value: list[str], attrs: Mapping | None = None) -> Any:
        """Only render selected choices when url is set."""
        if self.url:
            # Filter to only selected values
            selected_choices = {str(c) for c in value if c}
            all_choices = copy.copy(self.choices)

            # temporarily monkeypatch self.choices before calling super()
            try:
                self.choices = [c for c in self.choices if str(c[0]) in selected_choices]
                result = super().optgroups(name, value, attrs)
            finally:
                self.choices = all_choices
            return result

        return super().optgroups(name, value, attrs)


class ModelSelect2Multiple(QuerySetSelectMixin, Select2WidgetMixin, forms.SelectMultiple):
    """SelectMultiple widget for QuerySet choices and Select2."""
