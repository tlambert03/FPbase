from typing import TYPE_CHECKING

from django.db import models
from model_utils.models import TimeStampedModel

from proteins.models.fluorophore import Fluorophore
from proteins.models.mixins import Authorable, Product

if TYPE_CHECKING:
    from proteins.models import DyeState
    from references.models import Reference  # noqa: F401


class Dye(Authorable, TimeStampedModel, Product):  # TODO: rename to SmallMolecule
    """
    Represents a distinct organic fluorophore or chromophore.
    This is the 'Parent' entity for reactive derivatives.
    """

    # --- Identification ---
    name = models.CharField(max_length=255, db_index=True)
    slug = models.SlugField(unique=True)

    # Synonyms allow users to find "FITC" when searching "Fluorescein"
    # synonyms = ArrayField(models.CharField(max_length=255), blank=True, default=list)

    # --- Structural Status (The "Regret-Proof" Field) ---
    # Allows entry of proprietary dyes without forcing a fake structure.
    # STRUCTURE_STATUS_CHOICES = [
    #     ("DEFINED", "Defined Structure"),
    #     ("PROPRIETARY", "Proprietary / Unknown Structure"),
    # ]
    # structural_status = models.CharField(max_length=20, choices=STRUCTURE_STATUS_CHOICES, default="DEFINED")

    # --- Chemical Graph Data (Nullable for Proprietary Dyes) ---
    # We prioritize MolBlock for rendering, InChIKey for de-duplication.
    # canonical_smiles = models.TextField(blank=True)
    # inchi = models.TextField(blank=True)
    # inchikey = models.CharField(max_length=27, blank=True, db_index=True)
    # molblock = models.TextField(blank=True, help_text="V3000 Molfile for precise rendering")

    # --- Hierarchy & Ontology ---
    # Handles FITC (Parent) vs 5-FITC (Child) relationship
    # parent_mixture_id: int | None
    # parent_mixture = models.ForeignKey["Dye | None"](
    #     "self",
    #     on_delete=models.SET_NULL,
    #     null=True,
    #     blank=True,
    #     related_name="isomers",
    # )

    # Automated classification (e.g., "Rhodamine", "Cyanine", "BODIPY")
    # Populated via ClassyFire API or manual curation
    # chemical_class = models.CharField(max_length=100, blank=True, db_index=True)

    # --- Intrinsic Physics ---
    # Critical for fluorogenic dyes (JF dyes, SiR-tubulin)
    # Describes the Lactone-Zwitterion equilibrium constant.
    # equilibrium_constant_klz = models.FloatField(
    #     null=True,
    #     blank=True,
    #     help_text="Equilibrium constant between non-fluorescent lactone and fluorescent zwitterion.",
    # )

    if TYPE_CHECKING:
        states = models.manager.RelatedManager["DyeState"]()
        # isomers: models.QuerySet["Dye"]
        # derivatives: models.QuerySet[ReactiveDerivative]
        # collection_memberships: models.QuerySet

    class Meta:
        ...
        # Enforce uniqueness only on defined structures to allow multiple proprietary entries
        # constraints = [
        #     models.UniqueConstraint(
        #         fields=["inchikey"],
        #         name="unique_defined_molecule",
        #         condition=models.Q(structural_status="DEFINED"),
        #     )
        # ]

    def save(self, *args, **kwargs):
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
class DyeState(Fluorophore):
    """Represents a SmallMolecule in a specific environmental context.

    This holds the actual spectral data.
    """

    dye_id: int
    dye = models.ForeignKey["Dye"](Dye, on_delete=models.CASCADE, related_name="states")

    # --- Context ---
    # name = models.CharField(max_length=255, help_text="e.g., 'Bound to DNA' or 'In Methanol'")
    # solvent = models.CharField(max_length=100, default="PBS")
    # ph = models.FloatField(default=7.4)

    # --- Environmental Categorization ---
    # Helps the UI decide which spectrum to show for a specific query.
    # ENVIRONMENT_CHOICES = []
    # environment = models.CharField(max_length=20, choices=ENVIRONMENT_CHOICES, default="FREE")

    # --- Logic ---
    # is_reference = models.BooleanField(
    #     default=False, help_text="If True, this is the default state shown on the dye summary card."
    # )

    if TYPE_CHECKING:
        fluorophore_ptr: Fluorophore  # added by Django MTI

    def save(self, *args, **kwargs):
        self.entity_type = Fluorophore.EntityTypes.DYE
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


# This section handles the commercial reality of purchasing dyes.
#
# It separates "Chemist Tools" (Reactive Dyes) from "Biologist Tools" (Antibody Conjugates).


# class ReactiveDerivative(models.Model):
#     """A sold product derived from the SmallMolecule.

#     e.g., 'Janelia Fluor 549 NHS Ester' or 'JF549-HaloTag Ligand'
#     These are the products users buy to perform conjugation
#     (e.g., NHS esters, Maleimides, HaloTag Ligands).
#     """

#     core_dye_id: int
#     core_dye = models.ForeignKey["Dye"](
#         Dye,
#         on_delete=models.CASCADE,
#         related_name="derivatives",
#     )

#     # --- Chemistry ---
#     REACTIVE_GROUP_CHOICES = [
#         ("NHS_ESTER", "NHS Ester"),
#         ("HALO_TAG", "HaloTag Ligand"),
#         ("SNAP_TAG", "SNAP-Tag Ligand"),
#         ("CLIP_TAG", "CLIP-Tag Ligand"),
#         ("MALEIMIDE", "Maleimide"),
#         ("AZIDE", "Azide"),
#         ("ALKYNE", "Alkyne"),
#         ("BIOTIN", "Biotin"),
#         ("OTHER", "Other"),
#     ]
#     reactive_group = models.CharField(max_length=10, choices=REACTIVE_GROUP_CHOICES)

#     # Specific structure of the linker/handle (distinct from core dye)
#     full_smiles = models.TextField(blank=True, help_text="Structure of the complete reactive molecule")
#     molecular_weight = models.FloatField()

#     # --- Vendor Info ---
#     vendor = models.CharField(max_length=100)
#     catalog_number = models.CharField(max_length=100)

#     def __str__(self) -> str:
#         return f"{self.core_dye.common_name} - {self.reactive_group} ({self.vendor} {self.catalog_number})"


# Architectural Decision: We do not create a new DyeState or SmallMolecule for every
# antibody conjugate. Instead, we use a relational link.
#
# When a user views "Alexa 488 Anti-CD4", the system pulls:
#   1. Biology from the Antibody model.
#   2. Physics from the SmallMolecule model (specifically, the DyeState where
#      environment='PROTEIN_BOUND').


# class AntibodyConjugate(models.Model):
#     """
#     A virtual entity representing a commercial antibody-dye pairing.
#     Does NOT store spectral data; inherits it from the Fluorophore.
#     """

#     # The "Ingredients"
#     fluorophore_id: int
#     fluorophore = models.ForeignKey["Dye"](SmallMolecule, on_delete=models.PROTECT)
#     antibody_id: int
#     antibody = models.ForeignKey["Antibody"]("Antibody", on_delete=models.PROTECT)  # Define Antibody model elsewhere

#     # Product Details
#     vendor = models.CharField(max_length=100)
#     catalog_number = models.CharField(max_length=100)

#     # Heterogeneity
#     f_p_ratio = models.FloatField(
#         null=True, blank=True, help_text="Average Fluorophore/Protein ratio (Degree of Labeling)"
#     )

#     def get_display_spectrum(self):
#         """
#         Logic to fetch the correct spectrum.
#         Prioritizes a specific 'Protein Bound' state if available.
#         """
#         protein_state = self.fluorophore.states.filter(environment="PROTEIN_BOUND").first()
#         if protein_state:
#             return protein_state
#         # Fallback to reference state if no specific protein-bound data exists
#         return self.fluorophore.get_primary_spectrum()
