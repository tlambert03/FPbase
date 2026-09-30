from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import models
from django.db.models.signals import post_delete
from django.dispatch import receiver

from proteins.models.fluorescence_data import AbstractFluorescenceData
from proteins.models.fluorophore import FluorState

if TYPE_CHECKING:
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
        previous_state = None
        if rebuild_cache and not self._state.adding:
            # if the measurement is being moved, its previous state needs a rebuild too
            previous_state = FluorState.objects.filter(measurements=self.pk).first()
        super().save(*args, **kwargs)
        # Keep the parent cache in sync unless explicitly disabled
        if rebuild_cache:
            self.state.rebuild_attributes()
            if previous_state and previous_state.pk != self.state_id:
                previous_state.rebuild_attributes()


@receiver(post_delete, sender=FluorescenceMeasurement)
def _rebuild_state_after_delete(sender, instance: FluorescenceMeasurement, **kwargs) -> None:
    # a signal rather than a delete() override, so that queryset and cascading
    # deletes (which bypass Model.delete) also keep the parent cache in sync
    if state := FluorState.objects.filter(pk=instance.state_id).first():
        state.rebuild_attributes()
