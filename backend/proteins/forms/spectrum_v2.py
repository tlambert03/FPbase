"""Enhanced spectrum submission form with client-side processing and multi-spectrum support."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, TypedDict

from django import forms
from django.apps import apps
from django.db import transaction
from django.db.models import F
from django.utils.text import slugify

from proteins.extrest.entrez import is_valid_doi
from proteins.models import Dye, DyeState, FluorState, Spectrum, State
from references.models import Reference

if TYPE_CHECKING:
    from django.contrib.auth.models import User
    from django.db.models import QuerySet


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
    owner_slug: str | None  # Protein slug for autocomplete categories
    scale_factor: float | None
    ph: float | None
    solvent: str | None
    peak_wave: int | None


MAX_SPECTRA_PER_SUBMISSION = 20
MAX_DATA_POINTS_PER_SPECTRUM = 2000


def _validate_spectrum_json(raw: str | bytes) -> list[SpectrumJSONData]:
    if not raw or raw == "[]":
        raise forms.ValidationError("No spectrum data provided.")

    try:
        spectra = json.loads(raw)
    except json.JSONDecodeError as e:
        raise forms.ValidationError(f"Invalid JSON: {e}") from e

    if not isinstance(spectra, list) or len(spectra) == 0:
        raise forms.ValidationError("Expected a non-empty array of spectra.")

    if len(spectra) > MAX_SPECTRA_PER_SUBMISSION:
        raise forms.ValidationError(
            f"Too many spectra ({len(spectra)}). "
            f"Maximum {MAX_SPECTRA_PER_SUBMISSION} per submission."
        )

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

        if len(data) > MAX_DATA_POINTS_PER_SPECTRUM:
            raise forms.ValidationError(
                f"Spectrum {i + 1} has too many data points ({len(data)}). "
                f"Maximum {MAX_DATA_POINTS_PER_SPECTRUM}."
            )

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


def _protein_state(owner_slug: str | None) -> State | None:
    """The state that receives spectra submitted for the protein with `owner_slug`."""
    states = State.objects.select_related("protein").filter(protein__slug=owner_slug)
    # the protein's default state, falling back to its oldest state if none is set
    return states.order_by(F("default_for").desc(nulls_last=True), "id").first()


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
        help_text="Upload CSV or TSV file.",
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
        spectra = _validate_spectrum_json(raw)
        for spec in spectra:
            if existing := self._existing_spectra(spec).first():
                raise forms.ValidationError(
                    f"{spec['owner']} already has a spectrum of type "
                    f"'{existing.get_subtype_display()}'."
                )
        return spectra

    def _existing_spectra(self, spec: SpectrumJSONData) -> QuerySet[Spectrum]:
        """Approved spectra that the submitted `spec` would duplicate."""
        qs = Spectrum.objects.all()
        category, name = spec["category"], spec["owner"].strip()
        if category == Spectrum.PROTEIN:
            state = _protein_state(spec.get("owner_slug"))
            return qs.filter(owner_fluor=state, subtype=spec["subtype"]) if state else qs.none()
        if category == Spectrum.DYE:
            return qs.filter(
                owner_fluor__dyestate__dye__slug=slugify(name),
                owner_fluor__name=FluorState.DEFAULT_NAME,
                subtype=spec["subtype"],
            )
        # filters, cameras, and lights have a single spectrum, whatever its status
        owner_field = self.OWNER_LOOKUP[category][0]
        return Spectrum.objects.all_objects().filter(**{f"{owner_field}__name": name})

    def clean_primary_reference(self) -> str:
        """Validate that the DOI is resolvable if provided."""
        doi = self.cleaned_data.get("primary_reference", "").strip()
        if doi and not is_valid_doi(doi):
            raise forms.ValidationError(
                f"Could not find a reference for DOI: {doi}. Please check that it is correct."
            )

        return doi

    def clean(self):
        """Validate that at least one of source or primary_reference is provided."""
        cleaned_data = super().clean()
        source = cleaned_data.get("source", "").strip()
        reference = cleaned_data.get("primary_reference", "").strip()

        # Check if user attempted to provide a reference (even if it failed validation)
        # by looking at the raw data, not just cleaned_data
        attempted_reference = self.data.get("primary_reference", "").strip()

        if not source and not reference and not attempted_reference:
            raise forms.ValidationError(
                "Please provide at least one of Source or Primary Reference."
            )

        return cleaned_data

    def _get_or_create_owner(self, category: str, owner_name: str, owner_slug: str | None = None):
        """Get or create owner objects based on category.

        Args:
            category: The spectrum category (protein, dye, filter, etc.)
            owner_name: Display name of the owner
            owner_slug: For proteins, this is the Protein.slug from Select2 autocomplete

        Returns:
            Tuple of (owner_fluor, owner_filter, owner_camera, owner_light)
        """
        owner_fluor = owner_filter = owner_camera = owner_light = None

        if category == Spectrum.PROTEIN:
            # For proteins, owner_slug is the Protein.slug from Select2 autocomplete
            if not owner_slug:
                raise forms.ValidationError(
                    f"Protein '{owner_name}' must be selected from the autocomplete dropdown."
                )
            owner_fluor = _protein_state(owner_slug)
            if owner_fluor is None:
                raise forms.ValidationError(f"Protein not found: {owner_name}")

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
            owner_slug = spec_data.get("owner_slug")

            owner_fluor, owner_filter, owner_camera, owner_light = self._get_or_create_owner(
                category, owner_name, owner_slug
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
