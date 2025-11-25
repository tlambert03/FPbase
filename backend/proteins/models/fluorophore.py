from __future__ import annotations

from typing import TYPE_CHECKING, Final

from django.db import models
from django.utils.text import slugify

from proteins.models.fluorescence_data import AbstractFluorescenceData

if TYPE_CHECKING:
    from typing import Self

    from django.db.models import QuerySet
    from django.db.models.manager import RelatedManager

    from proteins.models import Dye, DyeState, FluorescenceMeasurement, OcFluorEff, Protein, State
    from proteins.models.spectrum import D3Dict, Spectrum


class FluorophoreManager[T: models.Model](models.Manager):
    _queryset_class: type[QuerySet[T]]

    def notdark(self):
        return self.filter(is_dark=False)

    def with_spectra(self):
        return self.get_queryset().filter(spectra__isnull=False).distinct()


# The Canonical Parent (The Summary)
class Fluorophore(AbstractFluorescenceData):
    """The database table for 'Things That Glow'.

    Polymorphic Fluorophore Parent.

    While fluorophores support multiple measurements of fluorescence data,
    `Fluorophore` also inherits `AbstractFluorescenceData`, and the values accessible
    on this instance serve as the canonical (i.e. "cached", "published", "composited")
    fluorescence properties for this entity.

    Contains the 'Accepted/Cached' values for generic querying.
    Acts as the materialized view of the 'best' measurements.
    """

    class EntityTypes(models.TextChoices):
        PROTEIN = ("p", "Protein")
        DYE = ("d", "Dye")

    # Identity
    DEFAULT_NAME: Final = "default"

    # State label (distinguishes states within same parent: "default", "red", "green")
    name = models.CharField(max_length=100, default=DEFAULT_NAME, db_index=True)

    # Cached parent info (denormalized for search performance)
    owner_name = models.CharField(
        max_length=255,
        db_index=True,
        blank=True,
        default="",
        help_text="Protein/Dye name (cached for searching)",
    )
    owner_slug = models.SlugField(
        max_length=200,
        blank=True,
        default="",
        help_text="Protein/Dye slug (cached for URLs)",
    )

    # Unique identifier (typically {owner_slug}-{state_name})
    slug = models.SlugField(max_length=200, unique=True)
    entity_type = models.CharField(max_length=2, choices=EntityTypes, db_index=True)

    # Lineage Tracking
    # Maps field names to Measurement IDs. e.g., {'ex_max': 102, 'qy': 105}
    source_map = models.JSONField(default=dict, blank=True)
    # Admin override: Per-field pinned measurement IDs that won't be auto-updated.
    # e.g., {'qy': 105} means qy always comes from measurement 105, ignoring priority rules.
    pinned_source_map = models.JSONField(default=dict, blank=True)

    # Managers
    objects: FluorophoreManager[Self] = FluorophoreManager()

    if TYPE_CHECKING:
        spectra: RelatedManager[Spectrum]
        measurements: RelatedManager[FluorescenceMeasurement]
        oc_effs: RelatedManager[OcFluorEff]

        # these are not *guaranteed* to exist, they come from Django MTI
        dyestate: DyeState
        state: State

    class Meta:
        indexes = [
            models.Index(fields=["ex_max"], name="fluorophore_ex_max_idx"),
            models.Index(fields=["em_max"], name="fluorophore_em_max_idx"),
            models.Index(fields=["owner_name"], name="fluorophore_owner_name_idx"),
            models.Index(fields=["entity_type", "is_dark"], name="fluorophore_type_dark_idx"),
        ]

    def __str__(self):
        return self.label

    def save(self, *args, **kwargs):
        # Auto-generate slug from owner_slug + state name if not set
        if not self.slug and self.owner_slug:
            self.slug = slugify(f"{self.owner_slug}-{self.name}")
        super().save(*args, **kwargs)

    @property
    def label(self) -> str:
        """Human-readable display name: 'EGFP' or 'mEos3.2 (red)'."""
        if not self.owner_name:
            return self.name
        if self.name == self.DEFAULT_NAME:
            return self.owner_name
        return f"{self.owner_name} ({self.name})"

    def as_subclass(self) -> Self:
        """Downcast to the specific subclass instance."""
        for subclass_name in ["dyestate", "state"]:
            if hasattr(self, subclass_name):
                return getattr(self, subclass_name)
        return self  # Fallback to parent if no child found

    def rebuild_attributes(self) -> None:
        """The Compositing Engine.

        Aggregates all measurements to determine the current canonical values.

        Priority order (highest to lowest):
        1. Pinned overrides - Admin has explicitly pinned a measurement for a field
        2. Trusted measurements - is_trusted=True on FluorescenceMeasurement
        3. Primary reference - Measurement from the owner's primary_reference
        4. Most recent - Fallback by date_measured (most recent first)
        """
        from proteins.models import FluorescenceMeasurement

        measurable_fields = AbstractFluorescenceData.get_measurable_fields()
        new_values: dict[str, object] = {}
        new_source_map: dict[str, int] = {}

        # Get primary reference ID for priority sorting
        primary_ref_id = self._get_primary_reference_id()

        # 1. Handle pinned fields first (admin overrides)
        pinned_fields: set[str] = set()
        for field, measurement_id in self.pinned_source_map.items():
            if field not in measurable_fields:
                continue
            try:
                measurement = FluorescenceMeasurement.objects.get(id=measurement_id, fluorophore=self)
                val = getattr(measurement, field)
                if val is not None:
                    new_values[field] = val
                    new_source_map[field] = measurement.id
                    pinned_fields.add(field)
            except FluorescenceMeasurement.DoesNotExist:
                # Pinned measurement no longer exists, skip it
                pass

        # 2. Fetch all measurements for non-pinned fields
        measurements = list(self.measurements.select_related("reference").all())

        # Sort by priority: trusted > primary_reference > most recent date
        def measurement_priority(m: FluorescenceMeasurement) -> tuple[int, int, str]:
            # Lower tuple = higher priority (sorted ascending)
            trusted_score = 0 if m.is_trusted else 1
            primary_score = 0 if primary_ref_id and m.reference_id == primary_ref_id else 1
            # Use date_measured for sorting, fallback to empty string (sorts last)
            date_str = m.date_measured.isoformat() if m.date_measured else ""
            # Negate by reversing string comparison for descending date order
            return (trusted_score, primary_score, date_str)

        # Sort: trusted first, then primary ref, then most recent date
        measurements.sort(key=measurement_priority)
        # Reverse date sorting within groups (we want most recent first)
        # Re-sort with proper date handling
        measurements.sort(
            key=lambda m: (
                0 if m.is_trusted else 1,
                0 if primary_ref_id and m.reference_id == primary_ref_id else 1,
                # Invert date for descending order (most recent first)
                -(m.date_measured.toordinal() if m.date_measured else 0),
            )
        )

        # 3. Waterfall Logic: Find the first non-null value for each non-pinned field
        for field in measurable_fields:
            if field in pinned_fields:
                continue  # Already handled above

            for m in measurements:
                val = getattr(m, field)
                if val is not None:
                    new_values[field] = val
                    new_source_map[field] = m.id
                    break
            else:
                # No measurement had a value for this field
                # Use field default for non-nullable fields (e.g., is_dark)
                field_obj = self._meta.get_field(field)
                if field_obj.has_default():
                    new_values[field] = field_obj.get_default()
                else:
                    new_values[field] = None

        # 4. Update cached values
        for key, val in new_values.items():
            setattr(self, key, val)

        self.source_map = new_source_map
        self.save()

    @property
    def fluor_name(self) -> str:
        if hasattr(self, "protein"):
            return self.protein.name
        return self.name

    @property
    def abs_spectrum(self) -> Spectrum | None:
        spect = [f for f in self.spectra.all() if f.subtype == "ab"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return None

    @property
    def ex_spectrum(self) -> Spectrum | None:
        spect = [f for f in self.spectra.all() if f.subtype == "ex"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return self.abs_spectrum

    @property
    def em_spectrum(self) -> Spectrum | None:
        spect = [f for f in self.spectra.all() if f.subtype == "em"]
        if len(spect) > 1:
            raise AssertionError(f"multiple em spectra found for {self}")
        if len(spect):
            return spect[0]
        return None

    @property
    def twop_spectrum(self) -> Spectrum | None:
        spect = [f for f in self.spectra.all() if f.subtype == "2p"]
        if len(spect) > 1:
            raise AssertionError("multiple 2p spectra found")
        if len(spect):
            return spect[0]
        return None

    @property
    def bright_rel_egfp(self) -> float | None:
        if self.brightness:
            return self.brightness / 0.336
        return None

    @property
    def stokes(self) -> float | None:
        try:
            return self.em_max - self.ex_max
        except TypeError:
            return None

    def has_spectra(self) -> bool:
        if any([self.ex_spectrum, self.em_spectrum]):
            return True
        return False

    def ex_band(self, height=0.7) -> tuple[float, float] | None:
        if (spect := self.ex_spectrum) is not None:
            return spect.width(height)
        return None

    def em_band(self, height=0.7) -> tuple[float, float] | None:
        if (spect := self.em_spectrum) is not None:
            return spect.width(height)
        return None

    def within_ex_band(self, value, height=0.7) -> bool:
        if band := self.ex_band(height):
            minRange, maxRange = band
            if minRange < value < maxRange:
                return True
        return False

    def within_em_band(self, value, height=0.7) -> bool:
        if band := self.em_band(height):
            minRange, maxRange = band
            if minRange < value < maxRange:
                return True
        return False

    def d3_dicts(self) -> list[D3Dict]:
        return [spect.d3dict() for spect in self.spectra.all()]

    def get_absolute_url(self) -> str | None:
        # return the absolute url for the protein or dye that owns this fluorophore
        if owner := self._owner():
            return owner.get_absolute_url()
        return None

    def _owner(self) -> Dye | Protein | None:
        if hasattr(self, "dyestate"):
            return self.dyestate.dye
        if hasattr(self, "state"):
            return self.state.protein
        return None

    def _get_primary_reference_id(self) -> int | None:
        """Get the primary reference ID for the owner entity (Protein or Dye)."""
        if owner := self._owner():
            return owner.primary_reference_id
        return None
