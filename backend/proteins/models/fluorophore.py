from typing import TYPE_CHECKING, Literal

from django.db import models

from proteins.models.fluorescence_data import AbstractFluorescenceData

if TYPE_CHECKING:
    from django.db.models.manager import RelatedManager

    from proteins.models import Spectrum


# The Canonical Parent (The Summary)
class Fluorophore(AbstractFluorescenceData):
    """The database table for 'Things That Glow'.

    Polymorphic Fluorophore Parent.

    While fluorophores support multiple measurements of fluorescence data,
    `Fluorophore` also inherits `AbstractFluorescenceData`, and the values accessible
    on this instance serve as the canonical (cached, "published", "composited")
    fluorescence properties for this entity.

    Contains the 'Accepted/Cached' values for generic querying.
    Acts as the materialized view of the 'best' measurements.
    """

    class EntityTypes(models.TextChoices):
        PROTEIN = ("protein", "Protein")
        DYE = ("dye", "Dye")

    # Identity (Hoisted for performance)
    label = models.CharField(max_length=255, db_index=True)
    slug = models.SlugField(unique=True)
    entity_type = models.CharField(max_length=10, choices=EntityTypes, db_index=True)

    # Lineage Tracking
    # Maps field names to Measurement IDs. e.g., {'ex_max': 102, 'qy': 105}
    source_map = models.JSONField(default=dict, blank=True)

    if TYPE_CHECKING:
        spectra = RelatedManager["Spectrum"]()

    class Meta:
        indexes = [
            models.Index(fields=["ex_max"]),
            models.Index(fields=["em_max"]),
        ]

    def __str__(self):
        return self.label

    def rebuild_attributes(self):
        """The Compositing Engine.

        Aggregates all measurements to determine the current canonical values.
        """
        # 1. Fetch all measurements, sorted by priority:
        #    Curator Trusted > Primary Reference > Most Recent Date
        measurements = self.measurements.select_related("reference").order_by(
            "-is_trusted",
            # Assuming you have a helper to check if ref is primary for the owner
            # This part requires custom logic depending on if self is Protein or Dye
            # '-is_primary_ref',
            "-date_measured",
        )

        measurable_fields = AbstractFluorescenceData.get_measurable_fields()
        new_values = {}
        new_source_map = {}

        # 2. Waterfall Logic: Find the first non-null value for each field
        for field in measurable_fields:
            found_val = None
            found_source_id = None

            for m in measurements:
                val = getattr(m, field)
                if val is not None:
                    found_val = val
                    found_source_id = m.id
                    break

            new_values[field] = found_val
            if found_source_id:
                new_source_map[field] = found_source_id

        # 3. Update Cache
        for key, val in new_values.items():
            setattr(self, key, val)

        self.source_map = new_source_map

        # 4. Refresh Label (Hoisting)
        if hasattr(self, "protein_state"):
            ps = self.protein_state
            self.label = f"{ps.protein.name} ({ps.name})"
        elif hasattr(self, "dye_state"):
            ds = self.dye_state
            self.label = f"{ds.dye.name} ({ds.name})"

        self.save()

    @property
    def fluor_name(self) -> str:
        if hasattr(self, "protein"):
            return self.protein.name
        return self.name

    @property
    def abs_spectrum(self) -> "Spectrum | None":
        spect = [f for f in self.spectra.all() if f.subtype == "ab"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return None

    @property
    def ex_spectrum(self) -> "Spectrum | None":
        spect = [f for f in self.spectra.all() if f.subtype == "ex"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return self.abs_spectrum

    @property
    def em_spectrum(self) -> "Spectrum | None":
        spect = [f for f in self.spectra.all() if f.subtype == "em"]
        if len(spect) > 1:
            raise AssertionError(f"multiple em spectra found for {self}")
        if len(spect):
            return spect[0]
        return self.abs_spectrum

    @property
    def twop_spectrum(self) -> "Spectrum | None":
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

    def ex_band(self, height=0.7) -> tuple[float, float] | Literal[False]:
        return self.ex_spectrum.width(height)

    def em_band(self, height=0.7) -> tuple[float, float] | Literal[False]:
        return self.em_spectrum.width(height)

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

    def d3_dicts(self) -> list[dict]:
        return [spect.d3dict() for spect in self.spectra.all()]
