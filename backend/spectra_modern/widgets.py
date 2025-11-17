"""Custom form widgets for the spectra_modern module."""

from __future__ import annotations

import json

from django import forms
from django.urls import reverse_lazy


class AutocompleteSelectWidget(forms.Select):
    """A select widget that renders as an autocomplete select using TomSelect.

    This widget renders a standard <select> with a single empty option,
    then relies on JavaScript (via the autocomplete-select.js module) to
    enhance it into a searchable autocomplete.

    The widget adds data-* attributes to configure the autocomplete behavior.

    Parameters
    ----------
    search_url : str, optional
        URL endpoint for searching. Can be a string or reverse_lazy() result.
    search_fields : list[str], optional
        List of field names to search on (e.g., ["text", "protein_name"]).
    placeholder : str, optional
        Placeholder text for the select, by default "Type to search..."
    min_query_length : int, optional
        Minimum query length before triggering search, by default 2
    attrs : dict, optional
        Additional HTML attributes for the select element
    choices : iterable, optional
        Initial choices (usually left empty for autocomplete)

    Examples
    --------
    ```py
    from spectra_modern.widgets import AutocompleteSelectWidget

    class MyForm(forms.Form):
        protein = forms.ModelChoiceField(
            queryset=State.objects.select_related("protein"),
            widget=AutocompleteSelectWidget(
                search_url=reverse_lazy("proteins:state-autocomplete"),
                search_fields=["text", "protein_name", "state_name"],
                placeholder="Type to search proteins...",
            ),
        )
    ```

    In the template, you can render the field normally:
        {{ form.protein }}

    Then initialize the autocomplete in JavaScript:
    ```js
        import { createAutocompleteSelectFromDataset } from './lib/autocomplete-select.js'
        document.querySelectorAll('[data-search-url]').forEach(el => {
            createAutocompleteSelectFromDataset(el)
        })
    ```
    """

    def __init__(
        self,
        search_url: str | None = None,
        search_fields: list[str] | None = None,
        placeholder: str = "Type to search...",
        min_query_length: int = 2,
        attrs: dict | None = None,
        choices=(),
    ):
        super().__init__(attrs=attrs, choices=choices)
        self.search_url = search_url
        self.search_fields = search_fields or ["text"]
        self.placeholder = placeholder
        self.min_query_length = min_query_length

    def get_context(self, name, value, attrs):
        """Build context with autocomplete configuration."""
        context = super().get_context(name, value, attrs)

        # Add data attributes for JavaScript configuration
        widget_attrs = context["widget"]["attrs"]
        if self.search_url:
            widget_attrs["data-search-url"] = self.search_url
        if self.search_fields:
            widget_attrs["data-search-fields"] = json.dumps(self.search_fields)
        if self.placeholder:
            widget_attrs["data-placeholder"] = self.placeholder
        if self.min_query_length != 2:  # Only set if different from default
            widget_attrs["data-min-query-length"] = self.min_query_length
        return context

    def optgroups(self, name, value, attrs=None):
        """
        Override to render only selected option(s) + empty option.

        This prevents rendering thousands of options in the initial HTML.
        The autocomplete will load options dynamically via AJAX.
        """
        # Get the selected value(s) if any
        if value:
            # For autocomplete, we only render the currently selected option
            # The rest will be loaded via AJAX when user searches
            selected_choices = [c for c in self.choices if str(c[0]) == str(value)]
        else:
            selected_choices = []

        # Create a minimal choice list: empty option + selected option(s)
        empty_label = getattr(self, "empty_label", None) or "---------"
        minimal_choices = [("", empty_label)]
        minimal_choices.extend(selected_choices)

        # Temporarily replace full choice list with minimal one
        original_choices = self.choices
        self.choices = minimal_choices
        try:
            return super().optgroups(name, value, attrs)
        finally:
            # Restore original choices
            self.choices = original_choices


class StateAutocompleteWidget(AutocompleteSelectWidget):
    """Convenience widget specifically for State (protein) autocomplete."""

    def __init__(self, attrs=None, choices=()):
        super().__init__(
            search_url=reverse_lazy("proteins:state-autocomplete"),
            search_fields=["text", "protein_name", "state_name"],
            placeholder="Type to search proteins...",
            attrs=attrs,
            choices=choices,
        )
