from typing import TYPE_CHECKING

from django.db import models
from django.utils.text import slugify
from model_utils.models import TimeStampedModel

from proteins.models.fluorophore import FluorState
from proteins.models.mixins import Authorable, Product

if TYPE_CHECKING:
    from proteins.models import DyeState
    from references.models import Reference


class Dye(Authorable, TimeStampedModel, Product):  # TODO: rename to SmallMolecule
    """
    Represents a distinct organic fluorophore or chromophore.
    This is the 'Parent' entity for reactive derivatives.
    """

    # --- Identification ---
    name = models.CharField(max_length=255, db_index=True)
    slug = models.SlugField(max_length=200, unique=True)

    default_state_id: int | None
    default_state: models.ForeignKey["DyeState | None"] = models.ForeignKey(
        "DyeState",
        related_name="default_for",
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )

    primary_reference_id: int | None
    primary_reference: models.ForeignKey["Reference | None"] = models.ForeignKey(
        "references.Reference",
        related_name="primary_dyes",
        verbose_name="Primary Reference",
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
        help_text="The publication that introduced the dye",
    )

    if TYPE_CHECKING:
        states = models.manager.RelatedManager["DyeState"]()

    class Meta:
        abstract = False

    def __str__(self) -> str:
        return self.name

    def save(self, *args, **kwargs):
        self.slug = slugify(self.name)  # Always regenerate, like Protein.save()
        super().save(*args, **kwargs)
        # Update cached owner fields on all states when dye name/slug changes
        # These fields are cached on Fluorophore for query performance
        if self.pk:
            self.states.update(owner_name=self.name, owner_slug=self.slug)

    def get_primary_spectrum(self):
        """Returns the 'Reference' DyeState (e.g. Protein-bound) for display cards."""
        return self.states.filter(is_reference=True).first()


# Instead of storing spectral data directly in the SmallMolecule, we link it here.
# This allows us to store "Alexa 488 in PBS" and "Alexa 488 in Ethanol" as valid,
# separate datasets.
class DyeState(FluorState):
    """Represents a SmallMolecule in a specific environmental context.

    This holds the actual spectral data.
    """

    dye_id: int
    dye: models.ForeignKey["Dye"] = models.ForeignKey(Dye, on_delete=models.CASCADE, related_name="states")

    if TYPE_CHECKING:
        fluorophore_ptr: FluorState  # added by Django MTI

    def save(self, *args, **kwargs):
        self.entity_type = FluorState.EntityTypes.DYE
        # Cache parent dye info for efficient searching
        if self.dye_id:
            self.owner_name = self.dye.name
            self.owner_slug = self.dye.slug
        super().save(*args, **kwargs)

    def get_absolute_url(self):
        # For now, return a URL based on the dye slug
        # TODO: implement proper dye detail view
        from django.urls import reverse

        return reverse("proteins:spectra") + f"?owner={self.dye.slug}"
