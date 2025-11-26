from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import models

if TYPE_CHECKING:
    from proteins.models import Protein


class SnapGenePlasmid(models.Model):
    """Represents a plasmid from SnapGene's database.

    To "refresh" the data from snapgene's website, run the management command:
        python manage.py sync_snapgene_plasmids
    """

    plasmid_id = models.CharField(max_length=100, unique=True, db_index=True)
    name = models.CharField(max_length=200)
    description = models.TextField(blank=True)
    author = models.CharField(max_length=200, blank=True)
    size = models.IntegerField(null=True, blank=True, help_text="Size in base pairs")
    topology = models.CharField(max_length=50, blank=True)

    if TYPE_CHECKING:
        proteins: models.QuerySet[Protein]

    class Meta:
        ordering = ["name"]
        verbose_name = "SnapGene Plasmid"
        verbose_name_plural = "SnapGene Plasmids"

    def __str__(self) -> str:
        return self.name

    @property
    def url(self) -> str:
        """Return the URL to this plasmid on SnapGene's website."""
        return f"https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids/{self.plasmid_id}"
