"""Clean Django forms without JS-specific dependencies."""

from __future__ import annotations

import ast

from django import forms
from django.apps import apps
from django.core.exceptions import ObjectDoesNotExist
from django.utils.safestring import mark_safe
from django.utils.text import slugify

from proteins.models import Fluorophore, Spectrum, State
from proteins.util.helpers import zip_wave_data
from proteins.util.importers import text_to_spectra
from proteins.validators import validate_spectrum


class SpectrumFormField(forms.CharField):
    """Custom field for spectrum data with validation."""

    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={"rows": 4, "placeholder": "[[wavelength, value], ...]"})

    def __init__(self, *args, **kwargs):
        if "help_text" not in kwargs:
            kwargs["help_text"] = (
                "List of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]. File data takes precedence."
            )
        super().__init__(*args, **kwargs)


class ModernSpectrumForm(forms.ModelForm):
    """
    Modern spectrum submission form without JavaScript dependencies.

    This form uses only standard Django widgets and relies on data attributes
    for progressive enhancement on the frontend. It does not require any specific
    JS libraries like Select2, autocomplete-light, or crispy-forms.
    """

    # Lookup mapping for spectrum categories to owner model names
    OWNER_LOOKUP = {
        Spectrum.DYE: ("owner_dye", "Dye"),
        Spectrum.PROTEIN: ("owner_state", "State"),
        Spectrum.FILTER: ("owner_filter", "Filter"),
        Spectrum.CAMERA: ("owner_camera", "Camera"),
        Spectrum.LIGHT: ("owner_light", "Light"),
    }

    # Valid spectrum subtypes for each category
    VALID_SUBTYPES = {
        "d": ["ex", "ab", "em", "2p"],  # Dye
        "p": ["ex", "ab", "em", "2p"],  # Protein
        "l": ["pd"],  # Light
        "f": ["bp", "bx", "bm", "sp", "lp", "bs"],  # Filter
        "c": ["qe"],  # Camera
        "": [],
    }

    # Protein state field (only shown when category is protein)
    owner_state = forms.ModelChoiceField(
        required=False,
        label="Protein",
        queryset=State.objects.select_related("protein"),
        widget=forms.Select(
            attrs={
                "data-autocomplete-url": "",  # Will be set in template
                "data-category": "protein",
            }
        ),
    )

    # Generic owner name field (for non-protein categories)
    owner = forms.CharField(
        max_length=100,
        label="Owner Name",
        required=False,
        help_text="Name of protein, dye, filter, etc...",
        widget=forms.TextInput(attrs={"data-validate-url": ""}),  # Will be set in template
    )

    # Spectrum data (manual entry)
    spectral_data = SpectrumFormField(required=False, label="Spectrum Data")

    # File upload
    file = forms.FileField(
        required=False,
        label="File Upload",
        help_text="2 column CSV/TSV file with wavelengths in first column and data in second column",
        widget=forms.FileInput(attrs={"accept": ".csv,.tsv,.txt"}),
    )

    # Confirmation checkbox
    confirmation = forms.BooleanField(
        required=True,
        label=mark_safe(
            "I understand that I am adding a spectrum to the <em>public</em> "
            "FPbase spectra database, and confirm that I have verified the validity of the data"
        ),
    )

    class Meta:
        model = Spectrum
        fields = (
            "category",
            "subtype",
            "file",
            "ph",
            "source",
            "solvent",
            "owner_state",
            "owner",
        )
        widgets = {
            "category": forms.Select(attrs={"data-valid-subtypes": ""}),  # Will be JSON in template
            "ph": forms.NumberInput(attrs={"step": "0.1", "min": "0", "max": "14"}),
        }
        # Note: spectral_data is a custom field that maps to the model's 'data' field

    def __init__(self, *args, user=None, **kwargs):
        """Initialize form with user context."""
        self.user = user
        super().__init__(*args, **kwargs)

        # Add Bootstrap classes to all fields
        for field in self.fields.values():
            if isinstance(field.widget, forms.CheckboxInput):
                field.widget.attrs["class"] = "form-check-input"
            else:
                current_class = field.widget.attrs.get("class", "")
                field.widget.attrs["class"] = f"{current_class} form-control".strip()

    def clean(self):
        """Validate that either file or data is provided."""
        cleaned_data = super().clean()

        # Check which data source was selected
        data_source = self.data.get("data_source", "file") if self.data else "file"

        if data_source == "manual":
            if not cleaned_data.get("spectral_data"):
                raise forms.ValidationError("Please enter valid spectrum data.")
        else:
            if not self.files.get("file"):
                raise forms.ValidationError("Please select a file to upload.")

        return cleaned_data

    def clean_file(self):
        """Process uploaded file and convert to spectrum data."""
        if self.files:
            filetext = ""
            x = []
            y = []
            try:
                for chunk in self.files["file"].chunks():
                    try:
                        filetext += chunk.decode("utf-8")
                    except AttributeError:
                        filetext += chunk
                x, y, _headers = text_to_spectra(filetext)
                if not len(y):
                    self.add_error("file", "Did not find a data column in the provided file")
                if not len(x):
                    self.add_error("file", "Could not parse wavelengths from first column")
            except Exception:
                self.add_error(
                    "file",
                    "Sorry, could not parse spectrum from this file. "
                    "Is it a two column CSV with (wavelength, spectrum)?",
                )
            if not self.errors:
                self.cleaned_data["spectral_data"] = zip_wave_data(x, y[0])
                # Also update the raw POST data dict
                self.data: dict = self.data.copy()
                self.data["spectral_data"] = self.cleaned_data["spectral_data"]

    def clean_owner_state(self):
        """Validate that protein doesn't already have this spectrum type."""
        owner_state = self.cleaned_data.get("owner_state")
        stype = self.cleaned_data.get("subtype")
        if self.cleaned_data.get("category") == Spectrum.PROTEIN:
            spectra = Spectrum.objects.all_objects().filter(owner_state=owner_state, subtype=stype)
            if spectra.exists():
                first = spectra.first()
                self.add_error(
                    "owner_state",
                    forms.ValidationError(
                        "%(owner)s already has a%(n)s %(stype)s spectrum %(status)s",
                        params={
                            "owner": owner_state,
                            "stype": first.get_subtype_display().lower(),
                            "n": "n" if stype != Spectrum.TWOP else "",
                            "status": " (pending)" if first.status == Spectrum.STATUS.pending else "",
                        },
                        code="owner_exists",
                    ),
                )
        return owner_state

    def clean_owner(self):
        """Validate that owner with same category and name doesn't already exist."""
        owner = self.cleaned_data.get("owner")
        cat = self.cleaned_data.get("category")
        stype = self.cleaned_data.get("subtype")

        if cat == Spectrum.PROTEIN:
            return owner

        try:
            mod = apps.get_model("proteins", self.OWNER_LOOKUP[cat][1])
            obj = mod.objects.get(slug=slugify(owner))
        except ObjectDoesNotExist:
            return owner
        except KeyError as e:
            if not cat:
                raise forms.ValidationError("Category not provided") from e
            else:
                raise forms.ValidationError("Category not recognized") from e
        else:
            # Object exists, check if it has this type of spectrum
            exists = False
            if isinstance(obj, Fluorophore):
                if obj.spectra.filter(subtype=stype).exists():
                    exists = True
                    stype = obj.spectra.filter(subtype=stype).first().get_subtype_display()
            elif hasattr(obj, "spectrum") and obj.spectrum:
                exists = True
                stype = obj.spectrum.get_subtype_display()

            if exists:
                self.add_error(
                    "owner",
                    forms.ValidationError(
                        "A %(model)s with the name %(name)s already has a %(stype)s spectrum",
                        params={
                            "model": self.OWNER_LOOKUP[cat][1].lower(),
                            "name": obj.name,
                            "stype": stype,
                        },
                        code="owner_exists",
                    ),
                )
            else:
                return owner

    def save(self, commit=True):
        """Save spectrum with appropriate owner."""
        # Map spectral_data field to model's data field
        if "spectral_data" in self.cleaned_data:
            data = self.cleaned_data["spectral_data"]
            # Convert string to list if needed
            if isinstance(data, str):
                data = ast.literal_eval(data)
            self.instance.data = data

        cat = self.cleaned_data.get("category")
        if cat != Spectrum.PROTEIN:
            owner_model = apps.get_model("proteins", self.OWNER_LOOKUP[cat][1])
            owner_name = self.cleaned_data.get("owner")
            # Only set created_by if user is authenticated
            defaults = {"created_by": self.user} if self.user.is_authenticated else {}
            owner_obj, created = owner_model.objects.get_or_create(name=owner_name, defaults=defaults)
            if not created and self.user.is_authenticated:
                owner_obj.updated_by = self.user
                owner_obj.save()
            setattr(self.instance, self.OWNER_LOOKUP[cat][0], owner_obj)
        # Only set created_by if user is authenticated
        if self.user.is_authenticated:
            self.instance.created_by = self.user
        return super().save(commit=commit)
