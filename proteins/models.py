# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.contrib.auth import get_user_model
from django.db.models import Avg
from django.utils.text import slugify
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator
from django.urls import reverse
from model_utils.models import StatusModel, TimeStampedModel
from model_utils import Choices
import uuid as uuid_lib
import json
import re

from references.models import Reference
from .helpers import fetch_ipg_sequence, get_color_group, wave_to_hex, mless, get_base_name
from .validators import protein_sequence_validator, validate_mutation
from .fields import SpectrumField
from reversion.models import Version

from Bio import Entrez
Entrez.email = "talley_lambert@hms.harvard.edu"

User = get_user_model()


class Authorable(models.Model):
    created_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_author')
    updated_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_modifier')

    class Meta:
        abstract = True


class Organism(Authorable, TimeStampedModel):
    """ A class for the parental organism (species) from which the protein has been engineered  """

    # Attributes
    id          = models.PositiveIntegerField(primary_key=True, verbose_name='Taxonomy ID', help_text="NCBI Taxonomy ID")  # genbank protein accession number
    scientific_name = models.CharField(max_length=128, blank=True)
    division    = models.CharField(max_length=128, blank=True)
    common_name = models.CharField(max_length=128, blank=True)
    species     = models.CharField(max_length=128, blank=True)
    genus       = models.CharField(max_length=128, blank=True)
    rank        = models.CharField(max_length=128, blank=True)

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = u'Organism'
        ordering = ['scientific_name']

    def save(self, *args, **kwargs):
        pubmed_record = Entrez.read(Entrez.esummary(db='taxonomy', id=self.id, retmode='xml'))
        self.scientific_name = pubmed_record[0]['ScientificName']
        self.division = pubmed_record[0]['Division']
        self.common_name = pubmed_record[0]['CommonName']
        self.species = pubmed_record[0]['Species']
        self.genus = pubmed_record[0]['Genus']
        self.rank = pubmed_record[0]['Rank']
        super(Organism, self).save(*args, **kwargs)


class Protein(Authorable, StatusModel, TimeStampedModel):
    """ Protein class to store individual proteins, each with a unique AA sequence and name  """

    STATUS = Choices('pending', 'approved')

    MONOMER = 'm'
    DIMER = 'd'
    TANDEM_DIMER = 'td'
    WEAK_DIMER = 'wd'
    TETRAMER = 't'
    AGG_CHOICES = (
        (MONOMER, 'Monomer'),
        (DIMER, 'Dimer'),
        (TANDEM_DIMER, 'Tandem dimer'),
        (WEAK_DIMER, 'Weak dimer'),
        (TETRAMER, 'Tetramer'),
    )

    BASIC = 'b'
    PHOTOACTIVATABLE = 'pa'
    PHOTOSWITCHABLE = 'ps'
    PHOTOCONVERTIBLE = 'pc'
    TIMER = 't'
    OTHER = 'o'
    SWITCHING_CHOICES = (
        (BASIC, 'Basic'),
        (PHOTOACTIVATABLE, 'Photoactivatable'),
        (PHOTOSWITCHABLE, 'Photoswitchable'),
        (PHOTOCONVERTIBLE, 'Photoconvertible'),
        (TIMER, 'Timer'),
        (OTHER, 'Other'),
    )

    # Attributes
    uuid        = models.UUIDField(default=uuid_lib.uuid4, editable=False, unique=True)  # for API
    name        = models.CharField(max_length=128, help_text="Name of the fluorescent protein", db_index=True)
    slug        = models.SlugField(max_length=64, unique=True, help_text="URL slug for the protein")  # for generating urls
    base_name   = models.CharField(max_length=128)  # easily searchable "family" name
    aliases     = ArrayField(models.CharField(max_length=200), blank=True, null=True)
    chromophore = models.CharField(max_length=5, null=True, blank=True)
    seq         = models.CharField(max_length=1024, unique=True, blank=True, null=True,
                    help_text="Amino acid sequence (IPG ID is preferred)",
                    validators=[protein_sequence_validator])  # consider adding Protein Sequence validator
    genbank     = models.CharField(max_length=12, null=True, blank=True, unique=True, verbose_name='Genbank Accession')
    uniprot     = models.CharField(max_length=12, null=True, blank=True, unique=True, verbose_name='UniProtKB Accession')
    ipg_id      = models.CharField(max_length=12, null=True, blank=True, unique=True,
                    verbose_name='IPG ID', help_text="Identical Protein Group ID at Pubmed")  # identical protein group uid
    mw          = models.FloatField(null=True, blank=True, help_text="Molecular Weight",)  # molecular weight
    agg         = models.CharField(max_length=2, choices=AGG_CHOICES, blank=True, help_text="Oligomerization tendency",)
    switch_type = models.CharField(max_length=2, choices=SWITCHING_CHOICES, blank=True,
                    verbose_name='Type', help_text="Photoswitching type (basic if none)")
    blurb       = models.CharField(max_length=512, blank=True, help_text="Brief descriptive blurb",)

    # Relations
    parent_organism = models.ForeignKey(Organism, related_name='proteins', verbose_name="Parental organism", blank=True, null=True, help_text="Organism from which the protein was engineered",)
    primary_reference = models.ForeignKey(Reference, related_name='primary_proteins', verbose_name="Primary Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Preferably the publication that introduced the protein",)  # usually, the original paper that published the protein
    references = models.ManyToManyField(Reference, related_name='proteins', blank=True)  # all papers that reference the protein
    FRET_partner = models.ManyToManyField('self', symmetrical=False, through='FRETpair', blank=True)
    default_state = models.ForeignKey('State', related_name='default_for', blank=True, null=True, on_delete=models.SET_NULL)

    __original_ipg_id = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # store IPG_ID so that we know if it changes
        self.__original_ipg_id = self.ipg_id

    @property
    def mless(self):
        return mless(self.name)

    @property
    def _base_name(self):
        '''return core name of protein, stripping prefixes like "m" or "Tag"'''
        return get_base_name(self.name)

    @property
    def versions(self):
        return Version.objects.get_for_object(self)

    def last_approved_version(self):
        if self.status == 'approved':
            return self
        try:
            return Version.objects.get_for_object(self) \
                                  .filter(serialized_data__contains='"status": "approved"') \
                                  .first()
        except Exception:
            return None

    @property
    def all_spectra(self):
        a = False
        if self.states.exists():
            a = []
            for s in self.states.all():
                if s.ex_spectra is not None:
                    a.append(s.ex_spectra.data)
                if s.em_spectra is not None:
                    a.append(s.em_spectra.data)
        return a

    @property
    def additional_references(self):
        return self.references.exclude(id=self.primary_reference_id).order_by('-year')

    @property
    def color(self):
        try:
            return get_color_group(self.default_state.ex_max, self.default_state.em_max)[0]
        except Exception:
            return ''

    # Methods
    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse("proteins:protein-detail", args=[self.slug])

    def has_default(self):
        return bool(self.default_state)

    def has_spectra(self):
        for state in self.states.all():
            if state.has_spectra():
                return True
        return False

    def spectra(self):
        spectra = []
        for state in self.states.all():
            spectra.append(state.nvd3ex)
            spectra.append(state.nvd3em)
        return json.dumps(spectra)

    def set_state_and_type(self):
        # FIXME: should allow control of default states in form
        # if only 1 state, make it the default state
        if not self.default_state or self.default_state.is_dark:
            if self.states.count() == 1 and not self.states.first().is_dark:
                self.default_state = self.states.first()
            # otherwise use farthest red non-dark state
            elif self.states.count() > 1:
                self.default_state = self.states.exclude(is_dark=True).order_by('-em_max').first()

        if self.states.count() == 1:
            self.switch_type = self.BASIC
        elif self.states.count() > 1:
            if not self.transitions.count():
                self.switch_type = self.OTHER
            if self.transitions.count() == 1:
                if self.states.filter(is_dark=True).count():
                    self.switch_type = self.PHOTOACTIVATABLE
                else:
                    self.switch_type = self.PHOTOCONVERTIBLE
            elif self.transitions.count() > 1:
                self.switch_type = self.PHOTOSWITCHABLE

    def clean(self):

        errors = {}
        # Don't allow basic switch_types to have more than one state.
#        if self.switch_type == 'b' and self.states.count() > 1:
#            errors.update({'switch_type': 'Basic (non photoconvertible) proteins cannot have more than one state.'})
        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):
        # Don't allow protein sequences to have non valid amino acid letters:
        if self.seq:
            self.seq = "".join(self.seq.split()).upper()  # remove whitespace

        # if the IPG ID has changed... refetch the sequence
        if self.ipg_id != self.__original_ipg_id:
            s = fetch_ipg_sequence(uid=self.ipg_id)
            self.seq = s[1] if s else None

        self.slug = slugify(self.name)
        self.base_name = self._base_name
        self.set_state_and_type()

        super().save(*args, **kwargs)
        self.__original_ipg_id = self.ipg_id

    # Meta
    class Meta:
        ordering = ['name']


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
    ext_coeff   = models.IntegerField(blank=True, null=True,
                    validators=[MinValueValidator(0), MaxValueValidator(300000)],
                    help_text="Extinction Coefficient")  # extinction coefficient
    qy          = models.FloatField(null=True, blank=True, help_text="Quantum Yield",
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
    protein     = models.ForeignKey(Protein, related_name="states", help_text="The protein to which this state belongs", on_delete=models.CASCADE)

    # Managers
    objects = StatesManager()

    @property
    def local_brightness(self):
        """ brightness relative to spectral neighbors.  1 = average """
        if not (self.em_max and self.brightness):
            return 1
        em_min = self.em_max - 15
        em_max = self.em_max + 15
        B = State.objects.exclude(id=self.id).filter(em_max__gt=em_min, em_max__lt=em_max).aggregate(Avg('brightness'))
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
        if self.ex_spectra is None and self.em_spectra is None:
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
                spec = self.ex_spectra.nvd3Format()
                peak = self.ex_max
                color = self.ex_spectra.color
            elif spectrum == 'em':
                spec = self.em_spectra.nvd3Format()
                peak = self.em_max
                color = self.em_spectra.color
            else:
                return
        except Exception:
            return None
        d = {
            "key": self.protein.name + " " + spectrum,
            "values": spec,
            "peak": peak,
            "type": spectrum,
            "color": color,
            "area": True,
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
        super(State, self).save(*args, **kwargs)

    class Meta:
        verbose_name = u'State'
        unique_together = (("protein", "ex_max", "em_max", "ext_coeff", "qy"),)


class BleachMeasurement(Authorable, TimeStampedModel):
    rate      = models.FloatField(verbose_name='Bleach Rate', help_text="Photobleaching rate",)  # bleaching half-life
    power     = models.FloatField(null=True, blank=True, verbose_name='Illumination Power', help_text="Illumination power (W/cm2)",)
    modality  = models.CharField(max_length=100, blank=True, verbose_name='Illumination Modality', help_text="Type of microscopy/illumination used for measurement",)
    reference = models.ForeignKey(Reference, related_name='bleach_measurements', verbose_name="Measurement Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Reference where the measurement was made",)  # usually, the original paper that published the protein
    state     = models.ForeignKey(State, related_name='bleach_measurements', verbose_name="Protein (state)", help_text="The protein (state) for which this measurement was observed", on_delete=models.CASCADE)

    def __str__(self):
        return "{}: {}{}".format(
            self.state,
            '{} s'.format(self.rate) if self.rate else '',
            'with'.format(self.modality) if self.modality else '')


class StateTransition(TimeStampedModel):
    trans_wave = models.PositiveSmallIntegerField(
            blank=True,
            null=True,
            verbose_name='Transition Wavelength',
            help_text="Wavelength required",
            validators=[MinValueValidator(300), MaxValueValidator(1000)]
            )
    protein = models.ForeignKey(Protein,
        related_name='transitions',
        verbose_name="Protein Transitioning",
        help_text="The protein that demonstrates this transition",
        on_delete=models.CASCADE)
    from_state = models.ForeignKey(State,
            related_name='transitions_from',
            verbose_name="From state",
            help_text="The initial state ",
            on_delete=models.CASCADE)
    to_state = models.ForeignKey(State,
            related_name='transitions_to',
            verbose_name="To state",
            help_text="The state after transition",
            on_delete=models.CASCADE)

    def clean(self):
        errors = {}
        if self.from_state.protein != self.protein:
            errors.update({'from_state': '"From" state must belong to protein {}'.format(protein.name)})
        if self.to_state.protein != self.protein:
            errors.update({'to_state': '"To" state must belong to protein {}'.format(protein.name)})
        if errors:
            raise ValidationError(errors)

    def __str__(self):
        return "Transition: {} {} -> {}".format(self.protein.name,
            self.from_state.name, self.to_state.name)


class FRETpair(Authorable, TimeStampedModel):
    # relational class for FRET pairs to hold attributes about the pair

    # Attributes
    radius   = models.FloatField(blank=True, null=True)

    # Relations
    donor    = models.ForeignKey(Protein, null=False, blank=False, verbose_name='donor', related_name='FK_FRETdonor_protein')
    acceptor = models.ForeignKey(Protein, null=False, blank=False, verbose_name='acceptor', related_name='FK_FRETacceptor_protein')

    pair_references = models.ManyToManyField(Reference, related_name='FK_FRETpair_reference', blank=True)  # any additional papers that reference the FRET pair

    @property
    def name(self):
        return self.donor.name + '-' + self.acceptor.name

    @property
    def spectral_overlap(self):
        accEx  = self.acceptor.default_state.ex_spectra
        accEC  = self.acceptor.default_state.ext_coeff
        donEm  = self.donor.default_state.em_spectra
        # donQY  = self.donor.default_state.qy
        donCum = sum(donEm.y)
        minAcc = accEx.min_wave
        maxAcc = accEx.max_wave
        minEm  = donEm.min_wave
        maxEm  = donEm.max_wave

        startingwave = int(max(minAcc, minEm))
        endingwave = int(min(maxAcc, maxEm))

        A = accEx.wave_value_pairs()
        D = donEm.wave_value_pairs()
        overlap = [(pow(wave, 4) * A[wave] * accEC * D[wave] / donCum) for wave in range(startingwave, endingwave + 1)]

        return sum(overlap)

    def forsterDist(self, n=1.4, k=2. / 3.):
        return .2108 * (pow((k) * (pow(n, -4) * self.spectral_overlap), (1. / 6.)))

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = 'FRET Pair'


class ProteinCollection(TimeStampedModel):
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=512, blank=True)
    proteins = models.ManyToManyField(Protein, related_name='collection_memberships')
    owner = models.ForeignKey(User, blank=True, null=True, related_name='collections', verbose_name='Protein Collection', on_delete=models.SET_NULL,)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse("proteins:collection-detail", args=[self.id])


class Mutation(models.Model):
    parent      = models.ForeignKey(Protein, related_name='proteins', verbose_name="Parent Protein")
    mutations   = ArrayField(models.CharField(max_length=5), validators=[validate_mutation])
    InvalidAlignment = Exception

    def child_seq(self):
        outseq = list(self.parent.seq)
        for mut in self.mutations:
            q = re.search(r'(?P<pre>\D+)(?P<pos>\d+)(?P<post>\D+)', mut)
            if q:
                pos = int(q.groupdict()['pos']) - 1
                pre = q.groupdict()['pre']
                if outseq[pos] == pre:
                    outseq[pos] = q.groupdict()['post']
                else:
                    raise self.InvalidAlignment('mutation letter {} at position {} \
                        does not agree with parent peptide {}'.format(pre, pos, outseq[pos]))
        return ''.join(outseq)

    def clean(self):
        for mut in self.mutations:
            q = re.search(r'(?P<pre>\D+)(?P<pos>\d+)', mut)
            if q:
                pos = int(q.groupdict()['pos']) - 1
                pre = q.groupdict()['pre']
                if not self.parent.seq[pos] == pre:
                    raise ValidationError('mutation {} at position {} does \
                        not agree with parent sequence {}'.format(pre, pos, self.parent.seq[pos]))

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)




