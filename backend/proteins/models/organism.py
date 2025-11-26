from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel

from proteins.extrest import entrez
from proteins.models.mixins import Authorable

if TYPE_CHECKING:
    from proteins.models import Protein


class Organism(Authorable, TimeStampedModel):
    """A class for the parental organism (species) from which the protein has been engineered"""

    # Attributes
    id = models.PositiveIntegerField(
        primary_key=True, verbose_name="Taxonomy ID", help_text="NCBI Taxonomy ID"
    )  # genbank protein accession number
    scientific_name = models.CharField(max_length=128, blank=True)
    division = models.CharField(max_length=128, blank=True)
    common_name = models.CharField(max_length=128, blank=True)
    species = models.CharField(max_length=128, blank=True)
    genus = models.CharField(max_length=128, blank=True)
    rank = models.CharField(max_length=128, blank=True)

    if TYPE_CHECKING:
        proteins: models.QuerySet[Protein]

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = "Organism"
        ordering = ["scientific_name"]

    def get_absolute_url(self):
        return reverse("proteins:organism-detail", args=[self.pk])

    def save(self, *args, **kwargs):
        if info := entrez.get_organism_info(self.id):
            self.__dict__.update(info)

        super().save(*args, **kwargs)

    def url(self):
        return self.get_absolute_url()
