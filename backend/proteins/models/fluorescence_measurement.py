from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import models

from proteins.models.fluorescence_data import AbstractFluorescenceData

if TYPE_CHECKING:
    from proteins.models import FluorState
    from references.models import Reference


# The "evidence"
class FluorescenceMeasurement(AbstractFluorescenceData):
    """Raw data points from a specific reference."""

    state_id: int
    state: models.ForeignKey[FluorState] = models.ForeignKey(
        "FluorState", related_name="measurements", on_delete=models.CASCADE
    )
    reference_id: int | None
    reference: models.ForeignKey[Reference | None] = models.ForeignKey(
        "references.Reference",
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )

    # Metadata specific to the act of measuring
    conditions = models.TextField(blank=True, help_text="pH, solvent, temp, etc.")

    def save(self, *args, rebuild_cache: bool = True, **kwargs) -> None:
        # Allow opt-out of rebuild during bulk operations to avoid N+1 queries
        super().save(*args, **kwargs)
        # Keep the parent cache in sync unless explicitly disabled
        if rebuild_cache:
            self.state.rebuild_attributes()

    def delete(self, *args, **kwargs) -> None:
        f = self.state
        super().delete(*args, **kwargs)
        # Always keep the parent cache in sync
        f.rebuild_attributes()
