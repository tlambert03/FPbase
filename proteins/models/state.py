# -*- coding: utf-8 -*-
from django.db import models
from django.db.models import Avg
from django.utils.text import slugify
from django.core.validators import MaxValueValidator, MinValueValidator
from model_utils.models import TimeStampedModel

from ..fields import SpectrumField
from .mixins import Authorable
from ..helpers import wave_to_hex


class StatesManager(models.Manager):
    def notdark(self):
        return self.filter(is_dark=False)


class State(Authorable, TimeStampedModel):
    """ A class for the states that a given protein can be in (including spectra and other state-dependent properties)  """

    # Attributes
    name        = models.CharField(max_length=64, default='default')  # required
    slug        = models.SlugField(max_length=128, unique=True, help_text="Unique slug for the state")  # calculated at save
    is_dark     = models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluorescence",)
    ex_max      = models.PositiveSmallIntegerField(blank=True, null=True,
                    validators=[MinValueValidator(300), MaxValueValidator(900)], db_index=True)
    em_max      = models.PositiveSmallIntegerField(blank=True, null=True,
                    validators=[MinValueValidator(300), MaxValueValidator(1000)], db_index=True)
    ex_spectra  = SpectrumField(blank=True, null=True, help_text='List of [[wavelength, value],...] pairs')  # excitation spectra (list of x,y coordinate pairs)
    em_spectra  = SpectrumField(blank=True, null=True, help_text='List of [[wavelength, value],...] pairs')  # emission spectra (list of x,y coordinate pairs)
    twop_ex_spectra = SpectrumField(blank=True, null=True, help_text='List of [[wavelength, value],...] pairs')  # 2 photon cross-section
    twop_ex_max = models.PositiveSmallIntegerField(blank=True, null=True, verbose_name='Peak 2P excitation',
                    validators=[MinValueValidator(700), MaxValueValidator(1600)], db_index=True)
    twop_peakGM = models.FloatField(null=True, blank=True, verbose_name='Peak 2P cross-section of S0->S1 (GM)',
                    validators=[MinValueValidator(0), MaxValueValidator(200)])
    twop_qy     = models.FloatField(null=True, blank=True, verbose_name="2P Quantum Yield",
                    validators=[MinValueValidator(0), MaxValueValidator(1)])  # quantum yield
    ext_coeff   = models.IntegerField(blank=True, null=True,
                    validators=[MinValueValidator(0), MaxValueValidator(300000)],
                    verbose_name="Extinction Coefficient")  # extinction coefficient
    qy          = models.FloatField(null=True, blank=True, verbose_name="Quantum Yield",
                    validators=[MinValueValidator(0), MaxValueValidator(1)])  # quantum yield
    brightness  = models.FloatField(null=True, blank=True, editable=False)
    pka         = models.FloatField(null=True, blank=True, verbose_name='pKa',
                    validators=[MinValueValidator(2), MaxValueValidator(12)])  # pKa acid dissociation constant
    maturation  = models.FloatField(null=True, blank=True, help_text="Maturation time (min)",  # maturation half-life in min
                    validators=[MinValueValidator(0), MaxValueValidator(1600)])
    lifetime    = models.FloatField(null=True, blank=True, help_text="Lifetime (ns)",
                    validators=[MinValueValidator(0), MaxValueValidator(20)])  # fluorescence lifetime in nanoseconds

    # Relations
    transitions = models.ManyToManyField('State', related_name='transition_state', verbose_name="State Transitions", blank=True, through='StateTransition')  # any additional papers that reference the protein
    protein     = models.ForeignKey('Protein', related_name="states", help_text="The protein to which this state belongs", on_delete=models.CASCADE)

    # Managers
    objects = StatesManager()

    @property
    def local_brightness(self):
        """ brightness relative to spectral neighbors.  1 = average """
        if not (self.em_max and self.brightness):
            return 1
        B = State.objects.exclude(id=self.id).filter(
                em_max__around=self.em_max).aggregate(Avg('brightness'))
        try:
            v = round(self.brightness / B['brightness__avg'], 4)
        except TypeError:
            v = 1
        return v

    @property
    def bright_rel_egfp(self):
        try:
            return round(float(self.brightness) / .336, 1)
        except TypeError:
            return None

    @property
    def stokes(self):
        try:
            return self.em_max - self.ex_max
        except TypeError:
            return None

    @property
    def nvd3ex(self):
        return self.nvd3dict('ex')

    @property
    def nvd3em(self):
        return self.nvd3dict('em')

    @property
    def nvd32p(self):
        return self.nvd3dict('2p')

    @property
    def emhex(self):
        if self.em_max and self.em_max > 0:
            return wave_to_hex(self.em_max)
        else:
            return '#000'

    @property
    def exhex(self):
        if self.ex_max and self.ex_max > 0:
            return wave_to_hex(self.ex_max)
        else:
            return '#000'

    # Methods
    def has_spectra(self):
        if self.ex_spectra is None and self.em_spectra is None and self.twop_ex_spectra is None:
            return False
        return True

    def EC_scaled_excitation(self):
        return [[n[0], n[1] * self.ext_coeff] for n in self.ex_spectra.data]

    def QY_scaled_emission(self):
        return [[n[0], n[1] * self.qy] for n in self.em_spectra.data]

    def ex_band(self, height=0.7):
        return self.ex_spectra.width(height)

    def em_band(self, height=0.7):
        return self.em_spectra.width(height)

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

    def nvd3dict(self, spectrum):
        try:
            if spectrum == 'ex':
                spec = self.ex_spectra
                peak = self.ex_max
                color = self.ex_spectra.color
            elif spectrum == 'em':
                spec = self.em_spectra
                peak = self.em_max
                color = self.em_spectra.color
            elif spectrum == '2p':
                spec = self.twop_ex_spectra
                peak = self.twop_ex_max
                color = wave_to_hex(self.twop_ex_max)
            else:
                return
        except Exception:
            return None
        d = {
            "key": "{}{} {}".format(self.protein.name,
                " " + self.name if not self.name == 'default' else "", spectrum),
            "values": spec.nvd3Format(),
            "peak": peak,
            "type": spectrum,
            "color": color,
            "area": True,
            "minwave": spec.min_wave,
            "maxwave": spec.max_wave,
        }
        return d

    def __str__(self):
        return "{} ({} state)".format(self.protein.name, self.name)

    def __repr__(self):
        return "<State: {}>".format(self.slug)

    def save(self, *args, **kwargs):
        self.slug = self.protein.slug + '_' + slugify(self.name)
        if self.qy and self.ext_coeff:
            self.brightness = float(round(self.ext_coeff * self.qy / 1000, 2))
        super().save(*args, **kwargs)

    class Meta:
        verbose_name = u'State'
        unique_together = (("protein", "ex_max", "em_max", "ext_coeff", "qy"),)
