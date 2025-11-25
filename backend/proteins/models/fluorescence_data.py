from __future__ import annotations

from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from model_utils.models import TimeStampedModel

from proteins.models.mixins import Authorable
from proteins.util.helpers import wave_to_hex


class AbstractFluorescenceData(Authorable, TimeStampedModel, models.Model):
    """Defines the physics schema.

    Used by both the Measurement (Input) and Fluorophore (Output/Cache).
    """

    ex_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(300), MaxValueValidator(900)],
        db_index=True,
        help_text="Excitation maximum (nm)",
    )
    em_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(300), MaxValueValidator(1000)],
        db_index=True,
        help_text="Emission maximum (nm)",
    )
    emhex = models.CharField(max_length=7, blank=True)
    exhex = models.CharField(max_length=7, blank=True)

    # core properties
    ext_coeff = models.IntegerField(
        blank=True,
        null=True,
        verbose_name="Extinction Coefficient (M-1 cm-1)",
        validators=[MinValueValidator(0), MaxValueValidator(300000)],
    )
    qy = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Quantum Yield",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )
    brightness = models.FloatField(null=True, blank=True, editable=False)
    lifetime = models.FloatField(
        null=True,
        blank=True,
        help_text="Lifetime (ns)",
        validators=[MinValueValidator(0), MaxValueValidator(20)],
    )
    pka = models.FloatField(
        null=True,
        blank=True,
        verbose_name="pKa",
        validators=[MinValueValidator(2), MaxValueValidator(12)],
    )

    # two photon properties
    twop_ex_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        verbose_name="Peak 2P excitation",
        validators=[MinValueValidator(700), MaxValueValidator(1600)],
        db_index=True,
    )
    twop_peak_gm = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Peak 2P cross-section of S0->S1 (GM)",
        validators=[MinValueValidator(0), MaxValueValidator(200)],
    )
    twop_qy = models.FloatField(
        null=True,
        blank=True,
        verbose_name="2P Quantum Yield",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )

    # extra

    is_dark = models.BooleanField(
        default=False,
        verbose_name="Dark State",
        help_text="This state does not fluoresce",
    )

    class Meta:
        abstract = True

    def save(self, *args, **kwargs):
        if self.qy and self.ext_coeff:
            self.brightness = float(round(self.ext_coeff * self.qy / 1000, 2))

        self.emhex = "#000" if self.is_dark else wave_to_hex(self.em_max)
        self.exhex = wave_to_hex(self.ex_max)
        super().save(*args, **kwargs)

    @classmethod
    def get_measurable_fields(cls):
        """Return only fluorescence property field names.

        Excludes metadata fields like id, created, modified, etc.
        Only returns fields that represent actual fluorescence measurements.
        """
        MEASURABLE = {
            "ex_max",
            "em_max",
            "ext_coeff",
            "qy",
            "brightness",
            "lifetime",
            "pka",
            "twop_ex_max",
            "twop_peak_gm",
            "twop_qy",
            "is_dark",
        }
        return [f.name for f in cls._meta.fields if f.name in MEASURABLE]
