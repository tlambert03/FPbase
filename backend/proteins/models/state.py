from django.contrib.contenttypes.fields import GenericRelation
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.db.models import Avg
from django.utils.text import slugify

from ..util.helpers import wave_to_hex
from .mixins import Product
from .spectrum import SpectrumOwner


class Fluorophore(SpectrumOwner):
    ex_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(300), MaxValueValidator(900)],
        db_index=True,
    )
    em_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        validators=[MinValueValidator(300), MaxValueValidator(1000)],
        db_index=True,
    )
    twop_ex_max = models.PositiveSmallIntegerField(
        blank=True,
        null=True,
        verbose_name="Peak 2P excitation",
        validators=[MinValueValidator(700), MaxValueValidator(1600)],
        db_index=True,
    )
    ext_coeff = models.IntegerField(
        blank=True,
        null=True,
        verbose_name="Extinction Coefficient",
        validators=[MinValueValidator(0), MaxValueValidator(300000)],
    )  # extinction coefficient
    twop_peakGM = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Peak 2P cross-section of S0->S1 (GM)",
        validators=[MinValueValidator(0), MaxValueValidator(200)],
    )
    qy = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Quantum Yield",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )  # quantum yield
    twop_qy = models.FloatField(
        null=True,
        blank=True,
        verbose_name="2P Quantum Yield",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )  # quantum yield
    brightness = models.FloatField(null=True, blank=True, editable=False)
    pka = models.FloatField(
        null=True,
        blank=True,
        verbose_name="pKa",
        validators=[MinValueValidator(2), MaxValueValidator(12)],
    )  # pKa acid dissociation constant
    lifetime = models.FloatField(
        null=True,
        blank=True,
        help_text="Lifetime (ns)",
        validators=[MinValueValidator(0), MaxValueValidator(20)],
    )  # fluorescence lifetime in nanoseconds
    emhex = models.CharField(max_length=7, blank=True)
    exhex = models.CharField(max_length=7, blank=True)
    is_dark = models.BooleanField(
        default=False,
        verbose_name="Dark State",
        help_text="This state does not fluorescence",
    )

    class Meta:
        abstract = True

    def save(self, *args, **kwargs):
        if self.qy and self.ext_coeff:
            self.brightness = float(round(self.ext_coeff * self.qy / 1000, 2))

        self.emhex = "#000" if self.is_dark else wave_to_hex(self.em_max)
        self.exhex = wave_to_hex(self.ex_max)

        super().save(*args, **kwargs)

    @property
    def fluor_name(self):
        if hasattr(self, "protein"):
            return self.protein.name
        return self.name

    @property
    def abs_spectrum(self):
        spect = [f for f in self.spectra.all() if f.subtype == "ab"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return None

    @property
    def ex_spectrum(self):
        spect = [f for f in self.spectra.all() if f.subtype == "ex"]
        if len(spect) > 1:
            raise AssertionError(f"multiple ex spectra found for {self}")
        if len(spect):
            return spect[0]
        return self.abs_spectrum

    @property
    def em_spectrum(self):
        spect = [f for f in self.spectra.all() if f.subtype == "em"]
        if len(spect) > 1:
            raise AssertionError(f"multiple em spectra found for {self}")
        if len(spect):
            return spect[0]
        return self.abs_spectrum

    @property
    def twop_spectrum(self):
        spect = [f for f in self.spectra.all() if f.subtype == "2p"]
        if len(spect) > 1:
            raise AssertionError("multiple 2p spectra found")
        if len(spect):
            return spect[0]
        return None

    @property
    def bright_rel_egfp(self):
        if self.brightness:
            return self.brightness / 0.336
        return None

    @property
    def stokes(self):
        try:
            return self.em_max - self.ex_max
        except TypeError:
            return None

    def has_spectra(self):
        if any([self.ex_spectrum, self.em_spectrum]):
            return True
        return False

    def ex_band(self, height=0.7):
        return self.ex_spectrum.width(height)

    def em_band(self, height=0.7):
        return self.em_spectrum.width(height)

    def within_ex_band(self, value, height=0.7):
        if self.has_spectra():
            minRange, maxRange = self.ex_band(height)
            if minRange < value < maxRange:
                return True
        return False

    def within_em_band(self, value, height=0.7):
        if self.has_spectra():
            minRange, maxRange = self.em_band(height)
            if minRange < value < maxRange:
                return True
        return False

    def d3_dicts(self):
        return [spect.d3dict() for spect in self.spectra.all()]


class FluorophoreManager(models.Manager):
    def notdark(self):
        return self.filter(is_dark=False)

    def with_spectra(self):
        return self.get_queryset().filter(spectra__isnull=False).distinct()


class Dye(Fluorophore, Product):
    objects = FluorophoreManager()
    oc_eff = GenericRelation("OcFluorEff", related_query_name="dye")


class State(Fluorophore):
    DEFAULT_NAME = "default"

    """ A class for the states that a given protein can be in
    (including spectra and other state-dependent properties)  """
    name = models.CharField(max_length=64, default=DEFAULT_NAME)  # required
    maturation = models.FloatField(
        null=True,
        blank=True,
        help_text="Maturation time (min)",  # maturation half-life in min
        validators=[MinValueValidator(0), MaxValueValidator(1600)],
    )
    # Relations
    transitions = models.ManyToManyField(
        "State",
        related_name="transition_state",
        verbose_name="State Transitions",
        blank=True,
        through="StateTransition",
    )  # any additional papers that reference the protein
    protein = models.ForeignKey(
        "Protein",
        related_name="states",
        help_text="The protein to which this state belongs",
        on_delete=models.CASCADE,
    )
    oc_eff = GenericRelation("OcFluorEff", related_query_name="state")

    # Managers
    objects = FluorophoreManager()

    class Meta:
        verbose_name = "State"
        unique_together = (("protein", "ex_max", "em_max", "ext_coeff", "qy"),)

    def __str__(self):
        if self.name in (self.DEFAULT_NAME, "default"):
            return str(self.protein.name)
        return f"{self.protein.name} ({self.name})"

    def get_absolute_url(self):
        return self.protein.get_absolute_url()

    def makeslug(self):
        return f"{self.protein.slug}_{slugify(self.name)}"

    @property
    def local_brightness(self):
        """brightness relative to spectral neighbors.  1 = average"""
        if not (self.em_max and self.brightness):
            return 1
        B = State.objects.exclude(id=self.id).filter(em_max__around=self.em_max).aggregate(Avg("brightness"))
        try:
            v = round(self.brightness / B["brightness__avg"], 4)
        except TypeError:
            v = 1
        return v
