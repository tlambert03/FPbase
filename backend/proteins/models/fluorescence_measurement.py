from typing import TYPE_CHECKING

from django.db import models

from proteins.models.fluorescence_data import AbstractFluorescenceData

if TYPE_CHECKING:
    from proteins.models import Fluorophore  # noqa: F401
    from references.models import Reference  # noqa: F401


# The "evidence"
class FluorescenceMeasurement(AbstractFluorescenceData):
    """Raw data points from a specific reference."""

    fluorophore_id: int
    fluorophore = models.ForeignKey["Fluorophore"](
        "Fluorophore", related_name="measurements", on_delete=models.CASCADE
    )
    reference_id: int
    reference = models.ForeignKey["Reference"]("references.Reference", on_delete=models.CASCADE)

    # Metadata specific to the act of measuring
    date_measured = models.DateField(null=True, blank=True)
    conditions = models.TextField(blank=True, help_text="pH, solvent, temp, etc.")

    # Curator Override
    is_trusted = models.BooleanField(default=False, help_text="If True, this measurement overrides others.")

    def save(self, *args, **kwargs) -> None:
        super().save(*args, **kwargs)
        # Always keep the parent cache in sync
        self.fluorophore.rebuild_attributes()

    def delete(self, *args, **kwargs) -> None:
        f = self.fluorophore
        super().delete(*args, **kwargs)
        # Always keep the parent cache in sync
        f.rebuild_attributes()
