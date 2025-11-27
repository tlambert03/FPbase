"""Enhanced spectrum submission form with client-side processing and multi-spectrum support."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Field, Layout, Submit
from dal import autocomplete
from django import forms
from django.apps import apps
from django.db import transaction
from django.utils.safestring import mark_safe
from django.utils.text import slugify

from proteins.models import Dye, DyeState, FluorState, Spectrum, State

if TYPE_CHECKING:
    from django.contrib.auth.models import User


class SpectrumFormV2(forms.Form):
    """Enhanced spectrum submission form supporting multi-column file uploads.

    This form handles client-side processing of spectrum data. The JavaScript frontend
    parses CSV/TSV files, allows column selection, normalizes data, and sends processed
    spectra as JSON.
    """

    # Lookup for non-protein, non-dye categories (filter/camera/light)
    lookup = {
        Spectrum.FILTER: ("owner_filter", "Filter"),
        Spectrum.CAMERA: ("owner_camera", "Camera"),
        Spectrum.LIGHT: ("owner_light", "Light"),
    }

    # Hidden field containing JSON array of processed spectra
    # Structure: [{ "data": [[wave, value]...], "subtype": "ex", "peak_wave": 488, ... }, ...]
    spectra_json = forms.CharField(
        widget=forms.HiddenInput(),
        required=True,
        error_messages={"required": "Please upload a file and select columns."},
    )

    # File upload field (for initial parsing by JavaScript)
    file = forms.FileField(
        required=False,  # Not required on POST since JS processes it client-side
        label="Spectrum File",
        help_text="Upload CSV or TSV file. Any column layout is accepted.",
    )

    # Shared metadata fields (apply to all spectra in batch)
    category = forms.ChoiceField(
        choices=Spectrum.CATEGORIES,
        label="Category",
        help_text="Type of spectrum data",
    )

    # For proteins: Select2 autocomplete
    owner_fluor = forms.ModelChoiceField(
        required=False,
        label=mark_safe('Protein<span class="asteriskField">*</span>'),
        queryset=State.objects.select_related("protein"),
        widget=autocomplete.ModelSelect2(
            url="proteins:state-autocomplete",
            attrs={"data-theme": "bootstrap-5", "data-width": "100%"},
        ),
    )

    # For non-protein categories (dye, filter, camera, light)
    owner = forms.CharField(
        max_length=100,
        label=mark_safe(
            '<span class="owner-type">Owner</span> Name<span class="asteriskField">*</span>'
        ),
        required=False,
        help_text="Name of dye, filter, camera, or light source",
    )

    # Optional metadata
    ph = forms.FloatField(
        required=False,
        label="pH",
        help_text="pH of the solution (for biological spectra)",
    )
    solvent = forms.CharField(
        max_length=100,
        required=False,
        label="Solvent",
        help_text="Solvent or buffer used",
    )
    source = forms.CharField(
        max_length=200,
        required=False,
        label="Source",
        help_text="Citation or source of the data",
    )

    # Confirmation checkbox
    confirmation = forms.BooleanField(
        required=True,
        label=mark_safe(
            "<span class='small'>I understand that I am adding spectrum data to the "
            "<em>public</em> FPbase spectra database, and confirm that I have verified "
            "the validity of the data</span>"
        ),
    )

    def __init__(self, *args, **kwargs):
        self.user: User | None = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.attrs = {"id": "spectrum-form-v2", "enctype": "multipart/form-data"}
        self.helper.add_input(Submit("submit", "Submit Spectra", css_class="btn-primary"))
        self.helper.layout = Layout(
            Div("category"),
            Div(
                Div("owner_fluor", css_class="col-sm-6 col-xs-12 protein-owner hidden"),
                Div("owner", css_class="col-sm-6 col-xs-12 non-protein-owner"),
                css_class="row",
            ),
            Div("file"),
            # Column picker and preview sections are rendered by JavaScript
            Div(css_id="column-picker-container", css_class="mb-3", style="display: none;"),
            Div(css_id="spectra-preview-container", css_class="mb-3", style="display: none;"),
            Div(
                Div("ph", css_class="col-md-6 col-sm-12"),
                Div("solvent", css_class="col-md-6 col-sm-12"),
                css_class="row bio-fields",
            ),
            Div("source"),
            Field("spectra_json", type="hidden"),
            Field("confirmation", css_class="custom-checkbox"),
        )

    def clean_spectra_json(self) -> list[dict]:
        """Parse and validate the JSON array of processed spectra."""
        raw = self.cleaned_data.get("spectra_json", "")
        if not raw:
            raise forms.ValidationError("No spectrum data provided.")

        try:
            spectra = json.loads(raw)
        except json.JSONDecodeError as e:
            raise forms.ValidationError(f"Invalid JSON: {e}") from e

        if not isinstance(spectra, list) or len(spectra) == 0:
            raise forms.ValidationError("Expected a non-empty array of spectra.")

        # Validate each spectrum
        for i, spec in enumerate(spectra):
            if not isinstance(spec, dict):
                raise forms.ValidationError(f"Spectrum {i + 1} is not a valid object.")

            # Required fields
            if "data" not in spec:
                raise forms.ValidationError(f"Spectrum {i + 1} is missing 'data' field.")

            data = spec["data"]
            if not isinstance(data, list) or len(data) < 2:
                raise forms.ValidationError(f"Spectrum {i + 1} must have at least 2 data points.")

            # Validate data structure: [[wavelength, value], ...]
            for j, point in enumerate(data):
                if not isinstance(point, list) or len(point) != 2:
                    raise forms.ValidationError(
                        f"Spectrum {i + 1}, point {j + 1}: must be [wavelength, value]."
                    )
                if not all(isinstance(v, (int, float)) for v in point):
                    raise forms.ValidationError(
                        f"Spectrum {i + 1}, point {j + 1}: values must be numbers."
                    )

            # Required subtype
            if "subtype" not in spec:
                raise forms.ValidationError(f"Spectrum {i + 1} is missing 'subtype'.")

            if spec["subtype"] not in dict(Spectrum.SUBTYPE_CHOICES):
                raise forms.ValidationError(
                    f"Spectrum {i + 1} has invalid subtype: {spec['subtype']}"
                )

        return spectra

    def clean(self):
        """Validate category-specific requirements."""
        cleaned_data = super().clean()
        category = cleaned_data.get("category")

        if category == Spectrum.PROTEIN:
            if not cleaned_data.get("owner_fluor"):
                self.add_error("owner_fluor", "Please select a protein state.")
        elif category in (Spectrum.DYE, Spectrum.FILTER, Spectrum.CAMERA, Spectrum.LIGHT):
            if not cleaned_data.get("owner"):
                self.add_error("owner", "Please enter the owner name.")

        return cleaned_data

    @transaction.atomic
    def save(self) -> list[Spectrum]:
        """Create Spectrum objects for each processed spectrum.

        Returns:
            List of created Spectrum objects.
        """
        spectra_data = self.cleaned_data["spectra_json"]
        category = self.cleaned_data["category"]
        created_spectra = []

        # Determine owner based on category
        owner_fluor = None
        owner_filter = None
        owner_camera = None
        owner_light = None

        if category == Spectrum.PROTEIN:
            owner_fluor = self.cleaned_data.get("owner_fluor")
        elif category == Spectrum.DYE:
            # Create or get Dye and DyeState
            owner_name = self.cleaned_data.get("owner")
            dye, created = Dye.objects.get_or_create(
                slug=slugify(owner_name),
                defaults={"name": owner_name, "created_by": self.user},
            )
            if not created and self.user:
                dye.updated_by = self.user
                dye.save()

            dye_state, _ = DyeState.objects.get_or_create(
                dye=dye,
                name=FluorState.DEFAULT_NAME,
                defaults={"created_by": self.user},
            )
            owner_fluor = dye_state
        elif category in self.lookup:
            # Filter, Camera, Light
            model_name = self.lookup[category][1]
            owner_model = apps.get_model("proteins", model_name)
            owner_name = self.cleaned_data.get("owner")
            owner_obj, created = owner_model.objects.get_or_create(
                name=owner_name,
                defaults={"created_by": self.user},
            )
            if not created and self.user:
                owner_obj.updated_by = self.user
                owner_obj.save()

            if category == Spectrum.FILTER:
                owner_filter = owner_obj
            elif category == Spectrum.CAMERA:
                owner_camera = owner_obj
            elif category == Spectrum.LIGHT:
                owner_light = owner_obj

        # Create each spectrum
        for spec_data in spectra_data:
            spectrum = Spectrum(
                category=category,
                subtype=spec_data["subtype"],
                owner_fluor=owner_fluor,
                owner_filter=owner_filter,
                owner_camera=owner_camera,
                owner_light=owner_light,
                ph=self.cleaned_data.get("ph"),
                solvent=self.cleaned_data.get("solvent") or "",
                source=self.cleaned_data.get("source") or "",
                created_by=self.user,
                # Set status based on user permissions
                status=Spectrum.STATUS.approved
                if self.user and self.user.is_staff
                else Spectrum.STATUS.pending,
            )

            # Set data using the property setter (handles normalization, etc.)
            spectrum.data = spec_data["data"]

            # Override computed values if provided by client
            if spec_data.get("peak_wave"):
                spectrum.peak_wave = spec_data["peak_wave"]

            if spec_data.get("scale_factor"):
                spectrum.scale_factor = spec_data["scale_factor"]

            spectrum.full_clean()
            spectrum.save()
            created_spectra.append(spectrum)

        return created_spectra
