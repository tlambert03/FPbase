from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from model_utils.models import TimeStampedModel

from references.models import Reference

from .mixins import Authorable


class OSERMeasurement(Authorable, TimeStampedModel):
    percent = models.FloatField(
        verbose_name="Percent Normal Cells",
        help_text="Percentage of 'normal' looking cells",
        blank=True,
        null=True,
        validators=[MinValueValidator(0), MaxValueValidator(100)],
    )
    percent_stddev = models.FloatField(
        verbose_name="StdDev",
        help_text="Standard deviation of percent normal cells (if applicable)",
        blank=True,
        null=True,
        validators=[MinValueValidator(0), MaxValueValidator(100)],
    )
    percent_ncells = models.IntegerField(
        blank=True,
        null=True,
        verbose_name="Number of cells for percent measurement",
        help_text="Number of cells analyzed in percent normal for this FP",
    )
    oserne = models.FloatField(
        verbose_name="OSER/NE ratio",
        help_text="Ratio of OSER to nuclear envelope (NE) fluorescence intensities",
        blank=True,
        null=True,
        validators=[MinValueValidator(0), MaxValueValidator(100)],
    )
    oserne_stddev = models.FloatField(
        verbose_name="OSER/NE StdDev",
        help_text="Standard deviation of OSER/NE ratio (if applicable)",
        blank=True,
        null=True,
        validators=[MinValueValidator(0), MaxValueValidator(100)],
    )
    oserne_ncells = models.IntegerField(
        blank=True,
        null=True,
        verbose_name="Number of cells for OSER/NE measurement",
        help_text="Number of cells analyzed in OSER/NE this FP",
    )
    celltype = models.CharField(
        max_length=64,
        blank=True,
        verbose_name="Cell Type",
        help_text="e.g. COS-7, HeLa",
    )
    temp = models.FloatField(null=True, blank=True, verbose_name="Temperature")

    reference = models.ForeignKey(
        Reference,
        related_name="oser_measurements",
        verbose_name="Measurement Reference",
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
        help_text="Reference where the measurement was made",
    )  # usually, the original paper that published the protein
    protein = models.ForeignKey(
        "Protein",
        related_name="oser_measurements",
        verbose_name="Protein",
        help_text="The protein on which this measurement was made",
        on_delete=models.CASCADE,
    )

    def __str__(self):
        return f"OSER: {self.protein} in {self.reference.citation}"

    class Meta:
        unique_together = (("protein", "reference"),)
