from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel
from Bio import Entrez
from .mixins import Authorable

Entrez.email = "talley_lambert@hms.harvard.edu"


class Organism(Authorable, TimeStampedModel):
    """ A class for the parental organism (species) from which the protein has been engineered  """

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

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = u"Organism"
        ordering = ["scientific_name"]

    def get_absolute_url(self):
        return reverse("proteins:organism-detail", args=[self.pk])

    def save(self, *args, **kwargs):
        pubmed_record = Entrez.read(
            Entrez.esummary(db="taxonomy", id=self.id, retmode="xml")
        )
        self.scientific_name = pubmed_record[0]["ScientificName"]
        self.division = pubmed_record[0]["Division"]
        self.common_name = pubmed_record[0]["CommonName"]
        self.species = pubmed_record[0]["Species"]
        self.genus = pubmed_record[0]["Genus"]
        self.rank = pubmed_record[0]["Rank"]
        super().save(*args, **kwargs)

    def url(self):
        return self.get_absolute_url()
