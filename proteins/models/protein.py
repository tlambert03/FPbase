# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.contrib.auth import get_user_model
from django.db.models import Q, Count
from django.utils.text import slugify
from django.core.exceptions import ValidationError
from django.urls import reverse
from django.conf import settings
from django.contrib.postgres.search import TrigramSimilarity
from model_utils.models import StatusModel, TimeStampedModel
from model_utils.managers import QueryManager
from model_utils import Choices

import uuid as uuid_lib
import json
import re
import io
import sys
import os
from subprocess import run, PIPE


from references.models import Reference
from .mixins import Authorable
from ..util.helpers import get_color_group, mless, get_base_name
# from .extrest.entrez import fetch_ipg_sequence
from ..validators import protein_sequence_validator, validate_uniprot
from .. import util
from ..util.sequence import FPSeq
from reversion.models import Version

from Bio import Entrez
Entrez.email = settings.ADMINS[0][1]

User = get_user_model()


def findname(name):
    queries = [
        {'name__iexact': name},
        {'aliases__icontains': name},
        # {'name__icontains': name},
        {'name__iexact': re.sub(r' \((Before|Planar|wild).*', '', name)},
        {'aliases__icontains': re.sub(r' \((Before|Planar|wild).*', '', name)},
        # {'name__icontains': re.sub(r' \((Before|Planar).*', '', name)},
        {'name__iexact': name.strip('1')},
        # {'name__icontains': name.strip('1')},
    ]
    for query in queries:
        try:
            return Protein.objects.get(**query)
        except Exception:
            pass
    return None


class ProteinQuerySet(models.QuerySet):

    def fasta(self):
        seqs = self.exclude(seq__isnull=True).values('slug', 'seq')
        return io.StringIO("\n".join([">{slug}\n{seq}".format(**s) for s in seqs]))

    def to_tree(self, output='clw'):
        fasta = self.fasta()
        if sys.platform == 'darwin':
            binary = 'bin/muscle_osx'
        else:
            binary = 'bin/muscle_nix'
        cmd = [binary]
        # faster
        cmd += ['-maxiters', '2', '-diags', '-quiet', '-' + output]
        # make tree
        cmd += ['-cluster', 'neighborjoining', '-tree2', 'tree.phy']
        result = run(cmd, input=fasta.read(), stdout=PIPE, encoding='ascii')
        with open('tree.phy', 'r') as handle:
            newick = handle.read().replace('\n', '')
        os.remove('tree.phy')
        return result.stdout, newick


class ProteinManager(models.Manager):

    def get_queryset(self):
        return ProteinQuerySet(self.model, using=self._db)

    def with_counts(self):
        from django.db import connection
        with connection.cursor() as cursor:
            cursor.execute("""
                SELECT p.id, p.name, p.slug, COUNT(*)
                FROM proteins_protein p, proteins_state s
                WHERE p.id = s.protein_id
                GROUP BY p.id, p.name, p.slug
                ORDER BY COUNT(*) DESC""")
            result_list = []
            for row in cursor.fetchall():
                p = self.model(id=row[0], name=row[1], slug=row[2])
                p.num_states = row[3]
                result_list.append(p)
        return result_list

    def annotated(self):
        return self.get_queryset().annotate(Count('states'), Count('transitions'))

    # def with_spectra(self):
    #     pids = [s.protein.id for s in
    #         State.objects.exclude(ex_spectra=None, em_spectra=None).distinct('protein')]
    #     return self.get_queryset().filter(id__in=pids)

    def with_spectra(self):
        return self.get_queryset().filter(states__spectra__isnull=False).distinct()

    def find_similar(self, name, similarity=0.2):
        return self.get_queryset() \
            .annotate(similarity=TrigramSimilarity('name', name)) \
            .filter(similarity__gt=similarity).order_by('-similarity')


class SequenceField(models.CharField):

    def __init__(self, *args, **kwargs):
        kwargs['max_length'] = 1024
        kwargs['validators'] = [protein_sequence_validator]
        super().__init__(*args, **kwargs)

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        del kwargs["max_length"]
        del kwargs["validators"]
        return name, path, args, kwargs

    def from_db_value(self, value, expression, connection):
        if not value:
            return None
        return FPSeq(value)

    def to_python(self, value):
        if isinstance(value, FPSeq):
            return value
        if not value:
            return None
        return FPSeq(value)

    def get_prep_value(self, value):
        if not value:
            return None
        return str(value)


class Protein(Authorable, StatusModel, TimeStampedModel):
    """ Protein class to store individual proteins, each with a unique AA sequence and name  """

    STATUS = Choices('pending', 'approved', 'hidden')

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
        (OTHER, 'Multistate'),
    )

    BILIVERDIN = 'bv'
    FLAVIN = 'fl'
    COFACTOR_CHOICES = (
        (BILIVERDIN, 'Biliverdin'),
        (FLAVIN, 'Flavin'),
    )

    # Attributes
    uuid        = models.UUIDField(default=uuid_lib.uuid4, editable=False, unique=True)  # for API
    name        = models.CharField(max_length=128, help_text="Name of the fluorescent protein", db_index=True)
    slug        = models.SlugField(max_length=64, unique=True, help_text="URL slug for the protein")  # for generating urls
    base_name   = models.CharField(max_length=128)  # easily searchable "family" name
    aliases     = ArrayField(models.CharField(max_length=200), blank=True, null=True)
    chromophore = models.CharField(max_length=5, null=True, blank=True)
    seq_validated = models.BooleanField(default=False, help_text='Sequence has been validated by a moderator')
    # seq must be nullable because of uniqueness contraints
    seq         = SequenceField(unique=True, blank=True, null=True, verbose_name='Sequence',
                                help_text="Amino acid sequence (IPG ID is preferred)")
    pdb         = ArrayField(models.CharField(max_length=4), blank=True, null=True, verbose_name='Protein DataBank ID')
    genbank     = models.CharField(max_length=12, null=True, blank=True, unique=True, verbose_name='Genbank Accession', help_text="NCBI Genbank Accession")
    uniprot     = models.CharField(max_length=10, null=True, blank=True, unique=True, verbose_name='UniProtKB Accession', validators=[validate_uniprot])
    ipg_id      = models.CharField(max_length=12, null=True, blank=True, unique=True,
                    verbose_name='IPG ID', help_text="Identical Protein Group ID at Pubmed")  # identical protein group uid
    mw          = models.FloatField(null=True, blank=True, help_text="Molecular Weight",)  # molecular weight
    agg         = models.CharField(max_length=2, choices=AGG_CHOICES, blank=True, help_text="Oligomerization tendency",)
    oser        = models.FloatField(null=True, blank=True, help_text="OSER score",)  # molecular weight
    switch_type = models.CharField(max_length=2, choices=SWITCHING_CHOICES, blank=True,
                    verbose_name='Type', help_text="Photoswitching type (basic if none)")
    blurb       = models.CharField(max_length=512, blank=True, help_text="Brief descriptive blurb",)
    cofactor    = models.CharField(max_length=2, choices=COFACTOR_CHOICES, blank=True, help_text='Required for fluorescence')
    # evo_parent  = models.ForeignKey('self', related_name='evo_children', verbose_name="Parental protein", on_delete=models.SET_NULL, blank=True, null=True, help_text="Protein from which the protein was evolved",)
    # evo_mutations = ArrayField(models.CharField(max_length=5), validators=[validate_mutation], blank=True, null=True)

    # Relations
    parent_organism = models.ForeignKey('Organism', related_name='proteins', verbose_name="Parental organism", on_delete=models.SET_NULL, blank=True, null=True, help_text="Organism from which the protein was engineered",)
    primary_reference = models.ForeignKey(Reference, related_name='primary_proteins', verbose_name="Primary Reference", blank=True, null=True, on_delete=models.SET_NULL, help_text="Preferably the publication that introduced the protein",)  # usually, the original paper that published the protein
    references = models.ManyToManyField(Reference, related_name='proteins', blank=True)  # all papers that reference the protein
    # FRET_partner = models.ManyToManyField('self', symmetrical=False, through='FRETpair', blank=True)
    default_state = models.ForeignKey('State', related_name='default_for', blank=True, null=True, on_delete=models.SET_NULL)

    __original_ipg_id = None

    # managers
    objects = ProteinManager()
    visible = QueryManager(~Q(status='hidden'))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # store IPG_ID so that we know if it changes
        self.__original_ipg_id = self.ipg_id

    @property
    def mless(self):
        return mless(self.name)

    @property
    def description(self):
        return util.long_blurb(self)

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

    def mutations_to(self, other, **kwargs):
        if isinstance(other, Protein):
            other = other.seq
        if not (self.seq and other):
            return None
        return self.seq.mutations_to(other, **kwargs)

    def mutations_from(self, other, **kwargs):
        if not (self.seq and other.seq):
            return None
        return other.seq.mutations_to(self.seq, **kwargs)

    def has_spectra(self):
        for state in self.states.all():
            if state.has_spectra():
                return True
        return False

    def has_bleach_measurements(self):
        return self.states.filter(bleach_measurements__isnull=False).exists()

    def d3_spectra(self):
        spectra = []
        for state in self.states.all():
            spectra.extend(state.d3_dicts())
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
        if self.pdb:
            self.pdb = list(set(self.pdb))
            for item in self.pdb:
                if Protein.objects.exclude(id=self.id).filter(pdb__contains=[item]).exists():
                    p = Protein.objects.filter(pdb__contains=[item]).first()
                    errors.update({'pdb': 'PDB ID {} is already in use by protein {}'.format(item, p.name)})

        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):

        # if the IPG ID has changed... refetch the sequence
        # if self.ipg_id != self.__original_ipg_id:
        #    s = fetch_ipg_sequence(uid=self.ipg_id)
        #    self.seq = s[1] if s else None

        self.slug = slugify(self.name)
        self.base_name = self._base_name
        self.set_state_and_type()

        super().save(*args, **kwargs)
        self.__original_ipg_id = self.ipg_id

    # Meta
    class Meta:
        ordering = ['name']
