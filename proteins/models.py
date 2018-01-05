# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.auth import get_user_model
from references.models import Reference
from django.utils.text import slugify
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db.models import Lookup
from django.urls import reverse
from model_utils.models import StatusModel, TimeStampedModel
from model_utils import Choices
import uuid as uuid_lib
import json
import re

from .helpers import fetch_ipg_sequence, get_color_group
from .validators import protein_sequence_validator

from Bio import Entrez
from Bio.Alphabet.IUPAC import protein as protein_alphabet
Entrez.email = "talley_lambert@hms.harvard.edu"

User = get_user_model()


class Around(Lookup):
    lookup_name = 'around'

    def as_sql(self, compiler, connection):
        lhs, lhs_params = self.process_lhs(compiler, connection)
        rhs, rhs_params = self.process_rhs(compiler, connection)
        params = lhs_params + rhs_params + lhs_params + rhs_params
        return '%s > %s - 16 AND %s < %s + 16' % (lhs, rhs, lhs, rhs), params


models.fields.PositiveSmallIntegerField.register_lookup(Around)


def wave_to_hex(wavelength, gamma=1):
    '''This converts a given wavelength into an approximate RGB value.
    The given wavelength is in nanometers.
    The range of wavelength is 380 nm through 750 nm.

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''

    wavelength = float(wavelength)
    if 520 <= wavelength:
        #pass
        wavelength += 40

    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 850:
        attenuation = 0.3 + 0.7 * (770 - wavelength) / (770 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    return '#%02x%02x%02x' % (int(R), int(G), int(B))


class Spectrum(object):
    """ Python class for spectra as a list of lists """

    def __init__(self, data=None):
        if data:
            if not isinstance(data, list):  # must be a list
                raise TypeError("Spectrum object must be of type List")
            if not all(isinstance(elem, list) for elem in data):  # must be list of lists
                raise TypeError("Spectrum object must be a list of lists")
            for elem in data:
                if not len(elem) == 2:
                    raise TypeError("All elements in Spectrum list must have two items")
                if not all(isinstance(n, (int, float)) for n in elem):
                    raise TypeError("All items in Septrum list elements must be numbers")
        self.data = data

    @property
    def x(self):
        self._x = []
        for i in self.data:
            self._x.append(i[0])
        return self._x

    @property
    def y(self):
        self._y = []
        for i in self.data:
            self._y.append(i[1])
        return self._y

    @property
    def color(self):
        return wave_to_hex(self.peak_wave)

    @property
    def peak_wave(self):
        return self.x[self.y.index(max(self.y))]

    @property
    def min_wave(self):
        return self.x[0]

    @property
    def max_wave(self):
        return self.x[-1]

    def __str__(self):
        return json.dumps(self.data)

    def width(self, height=0.5):
        try:
            upindex = next(x[0] for x in enumerate(self.y) if x[1] > height)
            downindex = len(self.y) - next(x[0] for x in enumerate(reversed(self.y)) if x[1] > height)
            return (self.x[upindex], self.x[downindex])
        except Exception:
            return False

    def change_x(self, value):
        if not isinstance(value, list):
            raise Exception("X values be a python list")
        if len(value) != len(self.data):
            raise Exception("Error: array length must match existing data")
        for i in range(len(value)):
            self.data[i][0] = value[i]

    def change_y(self, value):
        if not isinstance(value, list):
            raise Exception("Y values be a python list")
        if len(value) != len(self.data):
            raise Exception("Error: array length must match existing data")
        for i in range(len(value)):
            self.data[i][1] = value[i]

    def nvd3Format(self):
        output = []
        # arrayLength = len(self.data)
        for wave in range(350, int(self.min_wave)):
            output.append({'x': wave, 'y': 0})
        for elem in self.data:
            output.append({'x': elem[0], 'y': elem[1]})
        for wave in range(int(self.max_wave), 751):
            output.append({'x': wave, 'y': 0})
        return output

    def wave_value_pairs(self):
        output = {}
        # arrayLength = len(self.data)
        for elem in self.data:
            output[elem[0]] = elem[1]
        return output


class SpectrumField(models.TextField):
    description = "Stores a spectrum object"

    def __init__(self, *args, **kwargs):
        super(SpectrumField, self).__init__(*args, **kwargs)

    def from_db_value(self, value, expression, connection, context):
        if not value:
            return None
        return Spectrum(json.loads(value))

    def to_python(self, value):
        if isinstance(value, Spectrum):
            return value

        if not value:
            return None

        try:
            obj = json.loads(value)
            return Spectrum(obj)
        except Exception:
            raise ValidationError("Invalid input for a Spectrum instance")

    def get_prep_value(self, value):
        if value is None:
            return value
        return str(value)

    # def value_to_string(self, obj):
    #     value = self._get_val_from_obj(obj)
    #     return self.get_prep_value(value)


class Organism(TimeStampedModel):
    """ A class for the parental organism (species) from which the protein has been engineered  """

    # Attributes
    tax_id      = models.CharField(max_length=8, verbose_name='Taxonomy ID', help_text="NCBI Taxonomy ID (e.g. 6100 for Aequorea victora)",)  # genbank protein accession number
    scientific_name = models.CharField(max_length=128, blank=True)
    division    = models.CharField(max_length=128, blank=True)
    common_name = models.CharField(max_length=128, blank=True)
    species     = models.CharField(max_length=128, blank=True)
    genus       = models.CharField(max_length=128, blank=True)
    rank        = models.CharField(max_length=128, blank=True)

    # Relations
    created_by    = models.ForeignKey(User, related_name='organism_author', blank=True, null=True)  # the user who added the state
    updated_by  = models.ForeignKey(User, related_name='organism_modifier', blank=True, null=True)  # the user who last modified the state

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = u'Organism'
        ordering = ['scientific_name']

    def save(self, *args, **kwargs):
        pubmed_record = Entrez.read(Entrez.esummary(db='taxonomy', id=self.tax_id, retmode='xml'))
        self.scientific_name = pubmed_record[0]['ScientificName']
        self.division = pubmed_record[0]['Division']
        self.common_name = pubmed_record[0]['CommonName']
        self.species = pubmed_record[0]['Species']
        self.genus = pubmed_record[0]['Genus']
        self.rank = pubmed_record[0]['Rank']
        super(Organism, self).save(*args, **kwargs)


class Protein(TimeStampedModel):
    """ Protein class to store individual proteins, each with a unique AA sequence and name  """

#    STATUS = Choices('uncurated', 'curated', 'rejected')
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
    seq         = models.CharField(max_length=512, unique=True, blank=True, null=True,
                    help_text="Amino acid sequence (IPG ID is preferred)",
                    validators=[protein_sequence_validator])  # consider adding Protein Sequence validator
    gb_prot     = models.CharField(max_length=10, default='', blank=True, help_text="GenBank protein Accession number (e.g. AFR60231)",)  # genbank protein accession number
    gb_nuc      = models.CharField(max_length=10, default='', blank=True)  # genbank nucleotide accession number
    ipg_id      = models.CharField(max_length=12, null=True, blank=True, unique=True, verbose_name='IPG ID', help_text="Identical Protein Group ID at Pubmed")  # identical protein group uid
    mw          = models.FloatField(null=True, blank=True, help_text="Molecular Weight",)  # molecular weight
    agg         = models.CharField(max_length=2, choices=AGG_CHOICES, blank=True, help_text="Oligomerization tendency",)
    switch_type = models.CharField(max_length=2, choices=SWITCHING_CHOICES, blank=True, verbose_name='Type', help_text="Photoswitching type (basic if none)",)
    blurb       = models.CharField(max_length=512, blank=True, help_text="Brief descriptive blurb",)

    # Relations
    parent_organism = models.ForeignKey(Organism, related_name='proteins', verbose_name="Parental organism", blank=True, null=True, help_text="Organism from which the protein was engineered",)
    primary_reference = models.ForeignKey(Reference, related_name='primary_proteins', verbose_name="Primary Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Preferably the publication that introduced the protein",)  # usually, the original paper that published the protein
    references = models.ManyToManyField(Reference, related_name='proteins', verbose_name="References", blank=True)  # all papers that reference the protein
    FRET_partner = models.ManyToManyField('self', symmetrical=False, through='FRETpair', blank=True)
    created_by = models.ForeignKey(User, related_name='proteins_author', blank=True, null=True)  # the user who added the protein
    updated_by = models.ForeignKey(User, related_name='proteins_modifier', blank=True, null=True)
    default_state = models.ForeignKey('State', related_name='parent_protein', blank=True, null=True, on_delete=models.SET_NULL)

    __original_ipg_id = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__original_ipg_id = self.ipg_id

    # Manager
    # consider writing a manager to retrieve spectra or for custom queries

    @property
    def mless(self):
        name = self.name
        if re.search('^m[A-Z]', name):
            return name.lstrip('m')
        if name.startswith('monomeric'):
            name = name.lstrip('monomeric')
        if name.startswith('Monomeric'):
            name = name.lstrip('Monomeric')
        return name.lstrip(' ')

    @property
    def _base_name(self):
        '''return core name of protein, stripping prefixes like "m" or "Tag"'''
        name = self.name

        # remove PA/(Pa), PS, PC, from beginning
        if re.match('P[Aa]', name):
            name = name[2:]
        if re.match('P[Ss]', name):
            name = name[2:]
        if re.match('[Pp][Cc]', name):
            name = name[2:]
        if re.match('rs', name):
            name = name[2:]

        if re.match('LSS', name):
            name = name[3:].lstrip('-')

        # remove m (if next letter is caps) or monomeric
        if re.match('m[A-Z]', name):
            name = name[1:]
        if name.startswith('monomeric'):
            name = name.lstrip('monomeric')
        if name.startswith('Monomeric'):
            name = name.lstrip('Monomeric')

        # get rid of Td or td
        if re.match('[Tt][Dd][A-Z]', name):
            name = name[2:]
        if re.match('Tag', name):
            name = name[3:]

        # remove E at beginning (if second letter is caps)
        if re.match('E[A-Z]', name):
            name = name[1:]
        # remove S at beginning (if second letter is caps)
        if re.match('S[A-Z]', name):
            name = name[1:]

        # remove T- at beginning (if second letter is caps)
        if re.match('T-', name):
            name = name[2:]

        name = name.lstrip('-').lstrip(' ')

        return name

    @property
    def transitions(self):
        return [T for S in self.states.all() for T in S.transitions.all()]

    @property
    def count_states(self):
        return self.states.count()

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

    # Methods
    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse("proteins:protein-detail", args=[self.slug])

    # move these two methods to the state model
    def within_ex_band(self, value, height=0.7):
        if self.has_default():
            if self.default_state.has_spectra():
                minRange = self.default_state.ex_band(height)[0]
                maxRange = self.default_state.ex_band(height)[1]
                if minRange < value < maxRange:
                    return True
        return False

    def within_em_band(self, value, height=0.7):
        if self.has_default():
            if self.default_state.has_spectra():
                minRange = self.default_state.em_band(height)[0]
                maxRange = self.default_state.em_band(height)[1]
                if minRange < value < maxRange:
                    return True
        return False

    def has_default(self):
        if self.states.filter(default=True).count() > 0:
            return True
        else:
            return False

    def has_spectra(self):
        for state in self.states.all():
            if state.has_spectra():
                return True
        return False

    def color(self):
        try:
            return get_color_group(self.default_state.ex_max, self.default_state.em_max)[0]
        except Exception:
            return ''

    def spectra(self):
        return json.dumps([self.default_state.nvd3ex, self.default_state.nvd3em])

    def clean(self):
        # Don't allow protein sequences to have non valid amino acid letters:
        if self.seq:
            self.seq = "".join(self.seq.split()).upper()  # remove whitespace

        errors = {}
        # Don't allow basic switch_types to have more than one state.
#        if self.switch_type == 'b' and self.states.count() > 1:
#            errors.update({'switch_type': 'Basic (non photoconvertible) proteins cannot have more than one state.'})
        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):
        # if the IPG ID has changed... refetch the sequence
        if self.ipg_id != self.__original_ipg_id:
            s = fetch_ipg_sequence(uid=self.ipg_id)
            self.seq = s[1] if s else None

        self.slug = slugify(self.name)
        self.base_name = self._base_name

        # FIXME: should allow control of default states in form
        # if only 1 state, make it the default state
        try:
            self.default_state = self.states.get(default=True)
        except Exception:
            pass
        if not self.default_state:
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

        super().save(*args, **kwargs)
        self.__original_ipg_id = self.ipg_id

    # Meta
    class Meta:
        ordering = ['name']


class State(StatusModel, TimeStampedModel):
    """ A class for the states that a given protein can be in (including spectra and other state-dependent properties)  """
    STATUS = Choices('uncurated', 'curated', 'rejected')

    # Attributes
    name        = models.CharField(max_length=64, default='default')  # required
    slug        = models.SlugField(max_length=128, unique=True, help_text="Unique slug for the state")  # calculated at save
    is_dark     = models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluorescence",)
    ex_max      = models.PositiveSmallIntegerField(blank=True, null=True,
                    validators=[MinValueValidator(300), MaxValueValidator(900)], db_index=True)
    em_max      = models.PositiveSmallIntegerField(blank=True, null=True,
                    validators=[MinValueValidator(300), MaxValueValidator(1000)], db_index=True)
    ex_spectra  = SpectrumField(blank=True, null=True, help_text='Spectrum information as a list of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]')  # excitation spectra (list of x,y coordinate pairs)
    em_spectra  = SpectrumField(blank=True, null=True, help_text='Spectrum information as a list of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]')  # emission spectra (list of x,y coordinate pairs)
    ext_coeff   = models.IntegerField(blank=True, null=True,
                    validators=[MinValueValidator(0), MaxValueValidator(300000)],
                    help_text="Extinction Coefficient")  # extinction coefficient
    qy          = models.FloatField(null=True, blank=True, help_text="Quantum Yield",
                    validators=[MinValueValidator(0), MaxValueValidator(1)])  # quantum yield
    brightness  = models.FloatField(null=True, blank=True, editable=False)
    pka         = models.FloatField(null=True, blank=True, verbose_name='pKa',
                    validators=[MinValueValidator(2), MaxValueValidator(12)])  # pKa acid dissociation constant
    maturation  = models.FloatField(null=True, blank=True, help_text="Maturation time (min)",  # maturation half-life in min
                    validators=[MinValueValidator(0), MaxValueValidator(1400)])
    lifetime    = models.FloatField(null=True, blank=True, help_text="Lifetime (ns)",
                    validators=[MinValueValidator(0), MaxValueValidator(20)])  # fluorescence lifetime in nanoseconds

    # Relations
    transitions = models.ManyToManyField('State', related_name='transition_state', verbose_name="State Transitions", blank=True, through='StateTransition')  # any additional papers that reference the protein
    protein     = models.ForeignKey(Protein, related_name="states", help_text="The protein to which this state belongs", on_delete=models.CASCADE)
    created_by    = models.ForeignKey(User, related_name='state_author', blank=True, null=True)  # the user who added the state
    updated_by  = models.ForeignKey(User, related_name='state_modifier', blank=True, null=True)  # the user who last modified the state

    @property
    def local_brightness(self):
        """ brightness relative to spectral neighbors.  1 = average """
        from django.db.models import Avg
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
    def bleach(self):
        #TODO: this only gets the first bleaching measurement
        try:
            return self.bleach_measurement.first().rate
        except AttributeError:
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
            minRange = self.ex_band(height)[0]
            maxRange = self.ex_band(height)[1]
            if minRange < value < maxRange:
                return True
        return False

    def within_em_band(self, value, height=0.7):
        if self.has_spectra():
            minRange = self.em_band(height)[0]
            maxRange = self.em_band(height)[1]
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
        return self.slug

    def __repr__(self):
        return "<State: {}>".format(self.slug)

    class Meta:
        verbose_name = u'State'

    def save(self, *args, **kwargs):
        self.slug = self.protein.slug + '_' + slugify(self.name)
        if self.qy and self.ext_coeff:
            self.brightness = float(round(self.ext_coeff * self.qy / 1000, 2))
        super(State, self).save(*args, **kwargs)


class BleachMeasurement(TimeStampedModel):
    rate      = models.FloatField(verbose_name='Bleach Rate', help_text="Photobleaching rate",)  # bleaching half-life
    power     = models.FloatField(null=True, blank=True, verbose_name='Illumination Power', help_text="Illumination power (W/cm2)",)
    modality  = models.CharField(max_length=100, blank=True, verbose_name='Illumination Modality', help_text="Type of microscopy/illumination used for measurement",)
    reference = models.ForeignKey(Reference, related_name='bleach_measurement', verbose_name="Measurement Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Reference where the measurement was made",)  # usually, the original paper that published the protein
    state     = models.ForeignKey(State, related_name='bleach_measurement', verbose_name="Protein (state)", help_text="The protein (state) for which this measurement was observed", on_delete=models.CASCADE)


class StateTransition(TimeStampedModel):
    trans_wave = models.IntegerField(
            blank=True,
            null=True,
            verbose_name='Transition Wavelength',
            help_text="Wavelength of light that drives the protein into the 'To state:' (if applicable)"
            )
    protein = models.ForeignKey(Protein,
        related_name='transitions',
        verbose_name="Protein Transitioning",
        help_text="The protein that demonstrates this transition",
        on_delete=models.CASCADE)
    from_state = models.ForeignKey(State,
            related_name='transitions_from',
            verbose_name="From state",
            help_text="The initial state required for this transition to occur",
            on_delete=models.CASCADE)
    to_state = models.ForeignKey(State,
            related_name='transitions_to',
            blank=True,
            verbose_name="To state",
            help_text="The state to which this state transitions upon 'transition wavelength' illumination",
            on_delete=models.CASCADE)

    def __str__(self):
        return "<StateTransition: {} {}->{}>".format(self.protein.name,
            self.from_state, self.to_state)


# relational class for FRET pairs to hold attributes about the pair
class FRETpair(TimeStampedModel):

    # Attributes
    radius   = models.FloatField(blank=True, null=True)

    # Relations
    donor    = models.ForeignKey(Protein, null=False, blank=False, verbose_name='donor', related_name='FK_FRETdonor_protein')
    acceptor = models.ForeignKey(Protein, null=False, blank=False, verbose_name='acceptor', related_name='FK_FRETacceptor_protein')

    created_by    = models.ForeignKey(User, related_name="FRETpair_author", blank=True, null=True)  # the user who added the pair
    updated_by = models.ForeignKey(User, related_name='FRETpair_modifier', blank=True, null=True)  # the user who last modified the FRET pair data
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
        verbose_name = u'FRET Pair'