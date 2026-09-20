from __future__ import annotations

from typing import TYPE_CHECKING, Final, cast

from django.db import models
from django.db.models import Avg
from django.utils.text import slugify

from proteins.models.fluorescence_data import AbstractFluorescenceData

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Self

    from django.db.models import QuerySet
    from django.db.models.manager import RelatedManager

    from proteins.models import Dye, DyeState, FluorescenceMeasurement, OcFluorEff, Protein, State
    from proteins.models.spectrum import D3Dict, Spectrum


class FluorStateManager[T: models.Model](models.Manager):
    _queryset_class: type[QuerySet[T]]

    def notdark(self):
        return self.filter(is_dark=False)

    def with_spectra(self):
        return self.get_queryset().filter(spectra__isnull=False).distinct()


# The Canonical Parent (The Summary)
class FluorState(AbstractFluorescenceData):
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
    objects: FluorStateManager[Self] = FluorStateManager()

    if TYPE_CHECKING:
        spectra: RelatedManager[Spectrum]
        measurements: RelatedManager[FluorescenceMeasurement]
        oc_effs: RelatedManager[OcFluorEff]

        # these are not *guaranteed* to exist, they come from Django MTI
        dyestate: DyeState
        state: State

    class Meta:
        indexes = [
            models.Index(fields=["ex_max"], name="fluorstate_ex_max_idx"),
            models.Index(fields=["em_max"], name="fluorstate_em_max_idx"),
            models.Index(fields=["owner_name"], name="fluorstate_owner_name_idx"),
            models.Index(fields=["entity_type", "is_dark"], name="fluorstate_type_dark_idx"),
        ]

    def __str__(self):
        return self.label

    @classmethod
    def from_db(cls, db, field_names, values) -> Self:
        instance = super().from_db(db, field_names, values)
        instance._loaded_values = instance._measurable_values()
        return instance

    def _measurable_values(self) -> dict[str, object]:
        # (deferred fields are absent from __dict__; don't trigger a query for them)
        return {f: self.__dict__[f] for f in self._written_through_fields() if f in self.__dict__}

    @classmethod
    def _written_through_fields(cls) -> list[str]:
        # brightness is derived from ext_coeff and qy when a measurement is saved
        return [f for f in cls.get_measurable_fields() if f != "brightness"]

    def save(self, *args, write_through: bool = True, **kwargs):
        # Auto-generate slug from owner_slug + state name if not set
        if not self.slug and self.owner_slug:
            self.slug = slugify(f"{self.owner_slug}-{self.name}")
        super().save(*args, **kwargs)
        if write_through:
            # only fields that the caller changed: a state that is merely out of date
            # with respect to its measurements must not be written back over them
            # (a field that was deferred when loaded is only present if it was assigned)
            loaded = getattr(self, "_loaded_values", None)
            current = self._measurable_values()
            edited = [
                f
                for f, v in current.items()
                if loaded is None or f not in loaded or loaded[f] != v
            ]
            if (update_fields := kwargs.get("update_fields")) is not None:
                edited = [f for f in edited if f in update_fields]
            self.write_through(edited)
        self._loaded_values = self._measurable_values()

    def write_through(self, fields: Iterable[str] | None = None) -> list[str]:
        """Record values set directly on this state as a measurement, and rebuild.

        Measurements are the source of truth, so a value that lives only on the state
        would be lost by the next `rebuild_attributes`. Any of `fields` (default: all)
        that differ from what the measurements say are attributed to the owner's
        primary reference. Returns the names of the fields that were written.
        """
        if fields is None:
            fields = self._written_through_fields()
        if not (fields := list(fields)):
            return []
        composite, _ = self._composite()
        if not (fields := [f for f in fields if getattr(self, f) != composite[f]]):
            return []

        ref_id = self._get_primary_reference_id()
        measurement = self.measurements.filter(reference_id=ref_id).first()
        if measurement is None:
            author = self.updated_by or self.created_by
            measurement = self.measurements.model(
                state=self, reference_id=ref_id, created_by=author
            )
        for field in fields:
            setattr(measurement, field, getattr(self, field))
            # a direct edit supersedes a pin on that field
            self.pinned_source_map.pop(field, None)
        measurement.updated_by = self.updated_by
        measurement.save(rebuild_cache=False)
        self.rebuild_attributes()
        return fields

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
        """Set the canonical values (and `source_map`) from the measurements."""
        new_values, new_source_map = self._composite()
        for key, val in new_values.items():
            setattr(self, key, val)
        self.source_map = new_source_map
        self.save(write_through=False)

    def _composite(self) -> tuple[dict[str, object], dict[str, int]]:
        """The Compositing Engine.

        Aggregates all measurements to determine the current canonical values.
        Returns `(values, source_map)`.

        Priority order (highest to lowest):
        1. Pinned overrides - Admin has explicitly pinned a measurement for a field
        2. Primary reference - Measurement from the owner's primary_reference
        3. Others
        """
        measurable_fields = AbstractFluorescenceData.get_measurable_fields()
        new_values: dict[str, object] = {}
        new_source_map: dict[str, int] = {}
        primary_ref_id = self._get_primary_reference_id()

        # Fetch pinned measurements in one query
        pinned_source_map = cast("dict[str, int]", self.pinned_source_map)
        pinned_by_id = {
            m.id: m for m in self.measurements.filter(id__in=pinned_source_map.values())
        }

        # Handle pinned fields first (admin overrides)
        pinned_fields: set[str] = set()
        for field, mid in pinned_source_map.items():
            if (
                field in measurable_fields
                and (m := pinned_by_id.get(mid))
                and (val := getattr(m, field)) is not None
            ):
                new_values[field] = val
                new_source_map[field] = m.id
                pinned_fields.add(field)

        # Sort by primary_reference
        measurements = sorted(
            self.measurements.all(),
            # (an owner without a primary reference prefers unattributed measurements)
            key=lambda m: m.reference_id != primary_ref_id,
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

        return new_values, new_source_map

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

    @property
    def local_brightness(self) -> float:
        """Brightness relative to spectral neighbors. 1 = average."""
        if not (self.em_max and self.brightness):
            return 1
        avg = (
            FluorState.objects.exclude(id=self.id)
            .filter(em_max__around=self.em_max)
            .aggregate(Avg("brightness"))
        )
        try:
            return round(self.brightness / avg["brightness__avg"], 4)
        except TypeError:
            return 1

    def has_spectra(self) -> bool:
        return bool(any([self.ex_spectrum, self.em_spectrum]))

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
