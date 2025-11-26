from __future__ import annotations

from typing import TYPE_CHECKING, Final, cast

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
        4. others
        """
        measurable_fields = AbstractFluorescenceData.get_measurable_fields()
        new_values: dict[str, object] = {}
        new_source_map: dict[str, int] = {}
        primary_ref_id = self._get_primary_reference_id()

        # Fetch pinned measurements in one query
        pinned_source_map = cast("dict[str, int]", self.pinned_source_map)
        pinned_by_id = {m.id: m for m in self.measurements.filter(id__in=pinned_source_map.values())}

        # Handle pinned fields first (admin overrides)
        pinned_fields: set[str] = set()
        for field, mid in pinned_source_map.items():
            if field in measurable_fields and (m := pinned_by_id.get(mid)):
                if (val := getattr(m, field)) is not None:
                    new_values[field] = val
                    new_source_map[field] = m.id
                    pinned_fields.add(field)

        # Sort by priority: trusted > primary_reference
        measurements = sorted(
            self.measurements.all(),
            key=lambda m: (
                not m.is_trusted,
                not (primary_ref_id and m.reference_id == primary_ref_id),
            ),
        )

        # Waterfall: first non-null value for each non-pinned field
        for field in measurable_fields:
            if field in pinned_fields:
                continue
            for m in measurements:
                if (val := getattr(m, field)) is not None:
                    new_values[field] = val
                    new_source_map[field] = m.id
                    break
            else:
                field_obj = self._meta.get_field(field)
                new_values[field] = field_obj.get_default() if field_obj.has_default() else None

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
