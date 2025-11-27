"""Enhanced spectrum submission form with client-side processing and multi-spectrum support."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, TypedDict

from django import forms
from django.apps import apps
from django.db import transaction
from django.utils.text import slugify

from proteins.models import Dye, DyeState, FluorState, Spectrum, State
from references.models import Reference

if TYPE_CHECKING:
    from django.contrib.auth.models import User


class SpectrumJSONData(TypedDict):
    """Type definition for spectrum JSON data from frontend.

    Must match the SpectrumJSON typedef in form-controller.js.

    All fields are always present in the dict. Fields marked with | None
    can have null values when not applicable (e.g., ph/solvent for non-bio
    categories, scale_factor when not provided, peak_wave when not found).
    """

    # Required string/list fields (never None)
    data: list[list[float]]
    category: str
    owner: str
    subtype: str
    column_name: str

    # Always present but can be None
    scale_factor: float | None
    ph: float | None
    solvent: str | None
    peak_wave: int | None


def _validate_spectrum_json(raw: str | bytes) -> list[SpectrumJSONData]:
    if not raw or raw == "[]":
        raise forms.ValidationError("No spectrum data provided.")

    try:
        spectra = json.loads(raw)
    except json.JSONDecodeError as e:
        raise forms.ValidationError(f"Invalid JSON: {e}") from e

    if not isinstance(spectra, list) or len(spectra) == 0:
        raise forms.ValidationError("Expected a non-empty array of spectra.")

    valid_subtypes = dict(Spectrum.SUBTYPE_CHOICES)
    valid_categories = dict(Spectrum.CATEGORIES)

    for i, spec in enumerate(spectra):
        if not isinstance(spec, dict):
            raise forms.ValidationError(f"Spectrum {i + 1} is not a valid object.")

        # Validate data
        if "data" not in spec:
            raise forms.ValidationError(f"Spectrum {i + 1} is missing 'data' field.")

        data = spec["data"]
        if not isinstance(data, list) or len(data) < 2:
            raise forms.ValidationError(f"Spectrum {i + 1} must have at least 2 data points.")

        for j, point in enumerate(data):
            if not isinstance(point, list) or len(point) != 2:
                raise forms.ValidationError(
                    f"Spectrum {i + 1}, point {j + 1}: must be [wavelength, value]."
                )
            if not all(isinstance(v, (int, float)) for v in point):
                raise forms.ValidationError(
                    f"Spectrum {i + 1}, point {j + 1}: values must be numbers."
                )

        # Validate category
        if "category" not in spec or not spec["category"]:
            raise forms.ValidationError(f"Spectrum {i + 1} is missing category.")
        if spec["category"] not in valid_categories:
            raise forms.ValidationError(
                f"Spectrum {i + 1} has invalid category: {spec['category']}"
            )

        # Validate subtype
        if "subtype" not in spec or not spec["subtype"]:
            raise forms.ValidationError(f"Spectrum {i + 1} is missing subtype.")
        if spec["subtype"] not in valid_subtypes:
            raise forms.ValidationError(f"Spectrum {i + 1} has invalid subtype: {spec['subtype']}")

        # Validate owner
        if "owner" not in spec or not spec.get("owner", "").strip():
            raise forms.ValidationError(f"Spectrum {i + 1} is missing owner.")

    # Check for duplicate spectra within this submission
    # Use (category, owner, subtype) as the unique key
    seen = {}
    for i, spec in enumerate(spectra):
        key = (spec["category"], spec["owner"].strip().lower(), spec["subtype"])
        if key in seen:
            first_idx = seen[key]
            raise forms.ValidationError(
                f"Duplicate spectrum detected: Spectra {first_idx + 1} and {i + 1} have the same "
                f"owner ({spec['owner']}), category, and subtype ({spec['subtype']})."
            )
        seen[key] = i

    return spectra


class SpectrumFormV2(forms.Form):
    """Enhanced spectrum submission form supporting multi-spectrum file uploads.

    This form handles client-side processing of spectrum data. The JavaScript frontend
    parses CSV/TSV files, allows column selection, normalizes data, and sends processed
    spectra as JSON with per-spectrum metadata (category, owner, subtype, etc.).
    """

    # Lookup for non-protein, non-dye categories (filter/camera/light)
    OWNER_LOOKUP = {
        Spectrum.FILTER: ("owner_filter", "Filter"),
        Spectrum.CAMERA: ("owner_camera", "Camera"),
        Spectrum.LIGHT: ("owner_light", "Light"),
    }

    # Hidden field containing JSON array of processed spectra from JavaScript
    # Structure: [{ "data": [[wave, value]...], "category": "p", "owner": "EGFP",
    #              "subtype": "ex", "peak_wave": 488, ... }, ...]
    spectra_json = forms.CharField(
        widget=forms.HiddenInput(),
        required=True,
        error_messages={"required": "Please upload a file and configure your spectra."},
    )

    # File upload field (for initial parsing by JavaScript - not required on POST)
    file = forms.FileField(
        required=False,
        label="Spectrum File",
        help_text="Upload CSV or TSV file. Any column layout is accepted.",
    )

    # Shared source fields
    source = forms.CharField(
        max_length=200,
        required=False,
        label="Source",
        help_text="Citation or source of the data",
    )

    primary_reference = forms.CharField(
        max_length=200,
        required=False,
        label="Primary Reference (DOI)",
        help_text="Enter a valid DOI (e.g., 10.1234/example)",
    )

    # Confirmation checkbox
    confirmation = forms.BooleanField(
        required=True,
        label="I confirm the validity of this data",
    )

    def __init__(self, *args, **kwargs):
        self.user: User | None = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)

    def clean_spectra_json(self) -> list[SpectrumJSONData]:
        """Parse and validate the JSON array of processed spectra."""
        raw = self.cleaned_data.get("spectra_json", "")
        return _validate_spectrum_json(raw)

    def clean(self):
        """Validate that at least one of source or primary_reference is provided."""
        cleaned_data = super().clean()
        source = cleaned_data.get("source", "").strip()
        reference = cleaned_data.get("primary_reference", "").strip()

        if not source and not reference:
            raise forms.ValidationError(
                "Please provide at least one of Source or Primary Reference."
            )

        return cleaned_data

    def _get_or_create_owner(self, category: str, owner_name: str):
        """Get or create owner objects based on category.

        Returns:
            Tuple of (owner_fluor, owner_filter, owner_camera, owner_light)
        """
        owner_fluor = owner_filter = owner_camera = owner_light = None

        if category == Spectrum.PROTEIN:
            # Look up protein state by name
            try:
                owner_fluor = State.objects.select_related("protein").get(
                    protein__name__iexact=owner_name
                )
            except State.DoesNotExist:
                # Try by slug
                try:
                    owner_fluor = State.objects.select_related("protein").get(
                        protein__slug=slugify(owner_name)
                    )
                except State.DoesNotExist:
                    raise forms.ValidationError(f"Protein not found: {owner_name}") from None
            except State.MultipleObjectsReturned:
                # Get the default state
                owner_fluor = (
                    State.objects.select_related("protein")
                    .filter(protein__name__iexact=owner_name)
                    .first()
                )

        elif category == Spectrum.DYE:
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

        elif category in self.OWNER_LOOKUP:
            model_name = self.OWNER_LOOKUP[category][1]
            owner_model = apps.get_model("proteins", model_name)
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

        return owner_fluor, owner_filter, owner_camera, owner_light

    @transaction.atomic
    def save(self) -> list[Spectrum]:
        """Create Spectrum objects for each processed spectrum.

        Returns:
            List of created Spectrum objects.
        """
        spectra_data = self.cleaned_data["spectra_json"]
        source = self.cleaned_data.get("source", "")

        # Convert DOI string to Reference instance if provided
        reference_doi = self.cleaned_data.get("primary_reference", "").strip()
        reference = None
        if reference_doi:
            reference, _ = Reference.objects.get_or_create(doi=reference_doi)

        created_spectra = []

        for spec_data in spectra_data:
            category = spec_data["category"]
            owner_name = spec_data["owner"]

            owner_fluor, owner_filter, owner_camera, owner_light = self._get_or_create_owner(
                category, owner_name
            )

            spectrum = Spectrum(
                category=category,
                subtype=spec_data["subtype"],
                owner_fluor=owner_fluor,
                owner_filter=owner_filter,
                owner_camera=owner_camera,
                owner_light=owner_light,
                ph=spec_data.get("ph"),
                solvent=spec_data.get("solvent") or "",
                source=source,
                reference=reference,
                created_by=self.user,
                status=Spectrum.STATUS.approved
                if self.user and self.user.is_staff
                else Spectrum.STATUS.pending,
            )

            # Set data (handles normalization)
            spectrum.data = spec_data["data"]

            # Override computed values if provided
            if spec_data.get("peak_wave"):
                spectrum.peak_wave = spec_data["peak_wave"]
            if spec_data.get("scale_factor"):
                spectrum.scale_factor = spec_data["scale_factor"]

            spectrum.full_clean()
            spectrum.save()
            created_spectra.append(spectrum)

        return created_spectra
