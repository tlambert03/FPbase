# -*- coding: utf-8 -*-
from django.db import models
from django.conf import settings
from references.models import Reference
from django.template.defaultfilters import slugify
from django.core.exceptions import ValidationError
import json
import re

from Bio import Entrez
from Bio.Alphabet.IUPAC import protein as protein_alphabet
Entrez.email = "talley_lambert@hms.harvard.edu"

User = settings.AUTH_USER_MODEL

MONOMER = 'm'
DIMER = 'd'
TANDEM_DIMER = 'td'
WEAK_DIMER = 'wd'
TETRAMER = 't'
OLIGOMER_CHOICES = (
    (MONOMER, 'Monomer'),
    (DIMER, 'Dimer'),
    (TANDEM_DIMER, 'Tandem dimer'),
    (WEAK_DIMER, 'Weak dimer'),
    (TETRAMER, 'Tetramer'),
)


BASIC = 'b'
ACTIVATABLE = 'pa'
SWITCHABLE = 'ps'
CONVERTIBLE = 'pc'
TIMER = 't'
OTHER = 'o'
SWITCHING_CHOICES = (
    (BASIC, 'Basic'),
    (ACTIVATABLE, 'Photoactivatable'),
    (SWITCHABLE, 'Photoswitchable'),
    (CONVERTIBLE, 'Photoconvertible'),
    (TIMER, 'Timer'),
    (OTHER, 'Other'),
)


def wave_to_hex(wavelength, gamma=0.8):
    '''This converts a given wavelength into an approximate RGB value.
    The given wavelength is in nanometers.
    The range of wavelength is 380 nm through 750 nm.

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''

    wavelength = float(wavelength)

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
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
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
        arrayLength = len(self.data)
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
        except:
            raise ValidationError("Invalid input for a Spectrum instance")

    def get_prep_value(self, value):
        if value is None:
            return value
        return str(value)

    # def value_to_string(self, obj):
    #     value = self._get_val_from_obj(obj)
    #     return self.get_prep_value(value)


class Organism(models.Model):
    """ A class for the parental organism (species) from which the protein has been engineered  """

    # Attributes
    tax_id      = models.CharField(max_length=8, verbose_name='Taxonomy ID', help_text="Enter the NCBI Taxonomy ID (e.g. 6100 for Aequorea victora) and the info will be retrieved automatically",)  # genbank protein accession number
    scientific_name = models.CharField(max_length=128, blank=True)
    division    = models.CharField(max_length=128, blank=True)
    common_name = models.CharField(max_length=128, blank=True)
    species     = models.CharField(max_length=128, blank=True)
    genus       = models.CharField(max_length=128, blank=True)
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)

    # Relations
    added_by    = models.ForeignKey(User, related_name='organism_author', blank=True, null=True)  # the user who added the state
    updated_by  = models.ForeignKey(User, related_name='organism_modifiers', blank=True, null=True)  # the user who last modified the state

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = u'Organism'

    def save(self, *args, **kwargs):
        pubmed_record = Entrez.read(Entrez.esummary(db='taxonomy', id=self.tax_id, retmode='xml'))
        self.scientific_name = pubmed_record[0]['ScientificName']
        self.division = pubmed_record[0]['Division']
        self.common_name = pubmed_record[0]['CommonName']
        self.species = pubmed_record[0]['Species']
        self.genus = pubmed_record[0]['Genus']
        super(Organism, self).save(*args, **kwargs)


class ProteinManager(models.Manager):
    def get_queryset(self):
        return super(ProteinManager, self).get_queryset().filter(author='Roald Dahl')


class Protein(models.Model):
    """ Protein class to store individual proteins, each with a unique AA sequence and name  """

    # Attributes
    name        = models.CharField(max_length=128, help_text="Enter the name of the protein (required)", db_index=True)
    slug        = models.SlugField(max_length=64, unique=True, help_text="URL slug for the protein")  # for generating urls
    base_name   = models.CharField(max_length=128)  # easily searchable "family" name
    seq         = models.CharField(max_length=512, unique=True, blank=True, null=True, help_text="Amino acid sequence")  # consider adding Protein Sequence validator
    gb_prot     = models.CharField(max_length=10, null=True, blank=True, help_text="Enter the GenBank protein Accession number (e.g. AFR60231) and the sequence will be retrieved automatically",)  # genbank protein accession number
    gb_nuc      = models.CharField(max_length=10, null=True, blank=True)  # genbank nucleotide accession number
    ipg_id      = models.CharField(max_length=12, null=True, blank=True, unique=True, verbose_name='IPG ID', help_text="Identical Protein Group ID at Pubmed")  # identical protein group uid
    mw          = models.DecimalField(max_digits=5, decimal_places=2, null=True, blank=True, help_text="Molecular Weight",)  # molecular weight
    agg         = models.CharField(max_length=2, choices=OLIGOMER_CHOICES, blank=True, help_text="Oligomerization tendency",)
    switch_type = models.CharField(max_length=2, choices=SWITCHING_CHOICES, blank=True, verbose_name='Type', help_text="Photoswitching type (basic if none)",)
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)
    blurb       = models.CharField(max_length=512, blank=True, help_text="Brief descriptive blurb",)

    # Relations
    # default_state = models.ForeignKey('State', related_name='FK_defaultState_state', blank=True, null=True)  # default protein state
    parent_organism = models.ForeignKey(Organism, related_name='proteins', verbose_name="Parental organism", blank=True, null=True, help_text="Organism from which the protein was engineered",)
    primary_reference = models.ForeignKey(Reference, related_name='primary_proteins', verbose_name="Primary Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Preferably the publication that introduced the protein",)  # usually, the original paper that published the protein
    references = models.ManyToManyField(Reference, related_name='proteins', verbose_name="References", blank=True)  # all papers that reference the protein
    FRET_partner = models.ManyToManyField('self', symmetrical=False, through='FRETpair', blank=True)
    added_by = models.ForeignKey(User, related_name='proteins_author', blank=True, null=True)  # the user who added the protein
    updated_by = models.ForeignKey(User, related_name='proteins_modifiers', blank=True, null=True)
    default_state = models.ForeignKey('State', related_name='parent_protein', blank=True, null=True)

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

    # Methods
    def __str__(self):
        return self.name

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

    def clean(self):
        # Don't allow protein sequences to have non valid amino acid letters:
        errors = {}
        if self.seq:
            self.seq = "".join(self.seq.split()).upper()  # remove whitespace
            badletters = []
            for letter in self.seq:
                if letter not in protein_alphabet.letters:
                    badletters.append(letter)
            if len(badletters):
                badletters = set(badletters)
                errors.update({'seq': 'Invalid letter(s) found in amino acid sequence: {}'.format("".join(badletters))})

        # Don't allow basic switch_types to have more than one state.
        if self.switch_type == 'b' and self.states.count() > 1:
            errors.update({'switch_type': 'Basic (non photoconvertible) proteins cannot have more than one state.'})

        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):
        if self.gb_prot:
            pubmed_record = Entrez.read(Entrez.efetch(db='protein', id=self.gb_prot, retmode='xml'))
            if not self.seq:
                try:
                    self.seq = pubmed_record[0]['GBSeq_sequence'].upper()
                except Exception:
                    self.seq = None
        self.slug = slugify(self.name)
        self.base_name = self._base_name
        try:
            self.default_state = self.states.get(default=True)
        except Exception:
            pass
        self.full_clean()
        super(Protein, self).save(*args, **kwargs)

    # Meta
    class Meta:
        ordering = ['name']


class State(models.Model):
    """ A class for the states that a given protein can be in (including spectra and other state-dependent properties)  """

    # Attributes
    state_name  = models.CharField(max_length=128)  # required
    state_id    = models.CharField(max_length=128, unique=True)  # required
    default     = models.BooleanField(default=False, help_text="Check if this is the default (basal) state for the protein",)
    ex_max      = models.IntegerField(blank=True, null=True)
    em_max      = models.IntegerField(blank=True, null=True)
    ex_spectra  = SpectrumField(blank=True, null=True, help_text='Enter spectrum information as a list of [wavelength, value] pairs, e.g. [[300, 0.5],[301, 0.6] ... ]')  # excitation spectra (list of x,y coordinate pairs)
    em_spectra  = SpectrumField(blank=True, null=True, help_text='Enter spectrum information as a list of [wavelength, value] pairs, e.g. [[300, 0.5],[301, 0.6] ... ]')  # emission spectra (list of x,y coordinate pairs)
    ext_coeff   = models.IntegerField(blank=True, null=True, help_text="Extinction Coefficient",)  # extinction coefficient
    qy          = models.DecimalField(max_digits=4, decimal_places=3, null=True, blank=True, help_text="Quantum Yield")  # quantum yield
    pka         = models.DecimalField(max_digits=3, decimal_places=1, null=True, blank=True, verbose_name=u'pKa')  # pKa acid dissociation constant
#    bleach_wide = models.DecimalField(max_digits=5, decimal_places=1, null=True, blank=True, verbose_name='Bleach Widefield', help_text="Widefield photobleaching rate",)  # bleaching half-life for widefield microscopy
#    bleach_conf = models.DecimalField(max_digits=5, decimal_places=1, null=True, blank=True, verbose_name='Bleach Confocal', help_text="Confocal photobleaching rate",)  # bleaching half-life for confocal microscopy
    maturation  = models.DecimalField(max_digits=4, decimal_places=1, null=True, blank=True, help_text="Maturation time (seconds)")  # maturation half-life in minutes
    lifetime    = models.DecimalField(max_digits=3, decimal_places=2, null=True, blank=True, help_text="Fluorescence Lifetime (nanoseconds)",)  # fluorescence lifetime in nanoseconds
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)
    transitions = models.ManyToManyField('State', related_name='transition_state', verbose_name="State Transitions", blank=True, through='StateTransition')  # any additional papers that reference the protein
    # Relations
    protein     = models.ForeignKey(Protein, related_name="states", help_text="The protein to which this state belongs", on_delete=models.CASCADE)
    added_by    = models.ForeignKey(User, related_name='state_author', blank=True, null=True)  # the user who added the state
    updated_by  = models.ForeignKey(User, related_name='state_modifiers', blank=True, null=True)  # the user who last modified the state

    # Properties
    @property
    def brightness(self):
        try:
            return round(self.ext_coeff * self.qy / 1000, 2)
        except TypeError:
            return None

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
        return self.state_name

    def __repr__(self):
        return "<State: {}>".format(self.state_id)

    class Meta:
        verbose_name = u'State'

    def save(self, *args, **kwargs):
        self.state_id = self.protein.slug + "~" + slugify(self.state_name)
        super(State, self).save(*args, **kwargs)


class BleachMeasurement(models.Model):
    rate      = models.DecimalField(max_digits=6, decimal_places=1, verbose_name='Bleach Rate', help_text="Photobleaching rate",)  # bleaching half-life
    modality  = models.CharField(max_length=100, null=True, blank=True, verbose_name='Illumination Modality', help_text="Type of microscopy/illumination used for measurement",)
    reference = models.ForeignKey(Reference, related_name='bleach_measurement', verbose_name="Measurement Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Reference where the measurement was made",)  # usually, the original paper that published the protein
    state     = models.ForeignKey(State, related_name='bleach_measurement', verbose_name="Protein (state)", help_text="The protein (state) for which this measurement was observed", on_delete=models.CASCADE)


class StateTransition(models.Model):
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
class FRETpair(models.Model):

    # Attributes
    radius   = models.DecimalField(max_digits=4, decimal_places=2, blank=True, null=True)

    # Relations
    donor    = models.ForeignKey(Protein, null=False, blank=False, verbose_name='donor', related_name='FK_FRETdonor_protein')
    acceptor = models.ForeignKey(Protein, null=False, blank=False, verbose_name='acceptor', related_name='FK_FRETacceptor_protein')

    added_by    = models.ForeignKey(User, related_name="FRETpair_author", blank=True, null=True)  # the user who added the pair
    updated_by = models.ForeignKey(User, related_name='FRETpair_modifiers', blank=True, null=True)  # the user who last modified the FRET pair data
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)
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