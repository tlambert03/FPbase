import contextlib
import datetime
import io
import json
import os
import re
import sys
from random import choices
from subprocess import PIPE, run

import pytz
from Bio import Entrez
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.search import TrigramSimilarity
from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.db import models
from django.db.models import Count, F, Func, Q, Value
from django.urls import reverse
from django.utils.text import slugify
from model_utils import Choices
from model_utils.managers import QueryManager
from model_utils.models import StatusModel, TimeStampedModel
from reversion.models import Version

from favit.models import Favorite
from fpseq import FPSeq
from references.models import Reference

from .. import util
from ..util.helpers import get_base_name, get_color_group, mless, spectra_fig

# from .extrest.entrez import fetch_ipg_sequence
from ..validators import protein_sequence_validator, validate_uniprot
from .collection import ProteinCollection
from .mixins import Authorable
from .spectrum import Spectrum

Entrez.email = settings.ADMINS[0][1]

User = get_user_model()


def prot_uuid(k=5, opts="ABCDEFGHJKLMNOPQRSTUVWXYZ123456789"):
    i = "".join(choices(opts, k=k))
    try:
        Protein.objects.get(uuid=i)
    except Protein.DoesNotExist:
        return i
    else:
        return prot_uuid(k, opts)


def findname(name):
    queries = [
        {"name__iexact": name},
        {"aliases__icontains": name},
        # {'name__icontains': name},
        {"name__iexact": re.sub(r" \((Before|Planar|wild).*", "", name)},
        {"aliases__icontains": re.sub(r" \((Before|Planar|wild).*", "", name)},
        # {'name__icontains': re.sub(r' \((Before|Planar).*', '', name)},
        {"name__iexact": name.strip("1")},
        # {'name__icontains': name.strip('1')},
    ]
    for query in queries:
        with contextlib.suppress(Exception):
            return Protein.objects.get(**query)
    return None


class ProteinQuerySet(models.QuerySet):
    def fasta(self):
        seqs = list(self.exclude(seq__isnull=True).values("uuid", "name", "seq"))
        for s in seqs:
            s["name"] = s["name"].replace("\u03b1", "-alpha").replace("β", "-beta")
        return io.StringIO("\n".join([">{uuid} {name}\n{seq}".format(**s) for s in seqs]))

    def to_tree(self, output="clw"):
        fasta = self.fasta()
        binary = "bin/muscle_" + ("osx" if sys.platform == "darwin" else "nix")
        cmd = [binary]
        # faster
        cmd += ["-maxiters", "2", "-diags", "-quiet", f"-{output}"]
        # make tree
        cmd += ["-cluster", "neighborjoining", "-tree2", "tree.phy"]
        result = run(cmd, input=fasta.read(), stdout=PIPE, encoding="ascii")
        with open("tree.phy") as handle:
            newick = handle.read().replace("\n", "")
        os.remove("tree.phy")
        return result.stdout, newick


class ProteinManager(models.Manager):
    def deep_get(self, name):
        slug = slugify(name)
        with contextlib.suppress(Protein.DoesNotExist):
            return self.get(slug=slug)
        aliases_lower = Func(Func(F("aliases"), function="unnest"), function="LOWER")
        remove_space = Func(aliases_lower, Value(" "), Value("-"), function="replace")
        final = Func(remove_space, Value("."), Value(""), function="replace")
        D = dict(Protein.objects.annotate(aka=final).values_list("aka", "id"))
        if slug in D:
            return self.get(id=D[slug])
        else:
            raise Protein.DoesNotExist("Protein matching query does not exist.")

    def get_queryset(self):
        return ProteinQuerySet(self.model, using=self._db)

    def with_counts(self):
        from django.db import connection

        with connection.cursor() as cursor:
            cursor.execute(
                """
                SELECT p.id, p.name, p.slug, COUNT(*)
                FROM proteins_protein p, proteins_state s
                WHERE p.id = s.protein_id
                GROUP BY p.id, p.name, p.slug
                ORDER BY COUNT(*) DESC"""
            )
            result_list = []
            for row in cursor.fetchall():
                p = self.model(id=row[0], name=row[1], slug=row[2])
                p.num_states = row[3]
                result_list.append(p)
        return result_list

    def annotated(self):
        return self.get_queryset().annotate(Count("states"), Count("transitions"))

    # def with_spectra(self):
    #     pids = [s.protein.id for s in
    #         State.objects.exclude(ex_spectra=None, em_spectra=None).distinct('protein')]
    #     return self.get_queryset().filter(id__in=pids)

    def with_spectra(self, twoponly=False):
        qs = self.get_queryset().filter(states__spectra__isnull=False).distinct()
        if not twoponly:
            # hacky way to remove 2p only spectra
            qs = qs.annotate(stypes=Count("states__spectra__subtype")).filter(stypes__gt=1)
        return qs

    def find_similar(self, name, similarity=0.2):
        return (
            self.get_queryset()
            .annotate(similarity=TrigramSimilarity("name", name))
            .filter(similarity__gt=similarity)
            .order_by("-similarity")
        )


class SequenceField(models.CharField):
    def __init__(self, *args, **kwargs):
        kwargs["max_length"] = 1024
        kwargs["validators"] = [protein_sequence_validator]
        super().__init__(*args, **kwargs)

    def deconstruct(self):
        name, path, args, kwargs = super().deconstruct()
        del kwargs["max_length"]
        del kwargs["validators"]
        return name, path, args, kwargs

    def from_db_value(self, value, expression, connection):
        return FPSeq(value) if value else None

    def to_python(self, value):
        if isinstance(value, FPSeq):
            return value
        return FPSeq(value) if value else None

    def get_prep_value(self, value):
        return str(value) if value else None


# this is a hack to allow for reversions of proteins to work with Null chromophores
# this makes sure that a None value is converted to an empty string
class _NonNullChar(models.CharField):
    def to_python(self, value):
        return "" if value is None else super().to_python(value)


class Protein(Authorable, StatusModel, TimeStampedModel):
    """Protein class to store individual proteins, each with a unique AA sequence and name"""

    STATUS = Choices("pending", "approved", "hidden")

    MONOMER = "m"
    DIMER = "d"
    TANDEM_DIMER = "td"
    WEAK_DIMER = "wd"
    TETRAMER = "t"
    AGG_CHOICES = (
        (MONOMER, "Monomer"),
        (DIMER, "Dimer"),
        (TANDEM_DIMER, "Tandem dimer"),
        (WEAK_DIMER, "Weak dimer"),
        (TETRAMER, "Tetramer"),
    )

    BASIC = "b"
    PHOTOACTIVATABLE = "pa"
    PHOTOSWITCHABLE = "ps"
    PHOTOCONVERTIBLE = "pc"
    MULTIPHOTOCHROMIC = "mp"
    TIMER = "t"
    OTHER = "o"
    SWITCHING_CHOICES = (
        (BASIC, "Basic"),
        (PHOTOACTIVATABLE, "Photoactivatable"),
        (PHOTOSWITCHABLE, "Photoswitchable"),
        (PHOTOCONVERTIBLE, "Photoconvertible"),
        (MULTIPHOTOCHROMIC, "Multi-photochromic"),  # both convertible and switchable
        (OTHER, "Multistate"),
        (TIMER, "Timer"),
    )

    BILIRUBIN = "br"
    BILIVERDIN = "bv"
    FLAVIN = "fl"
    PHYCOCYANOBILIN = "pc"
    RIBITYL_LUMAZINE = "rl"
    COFACTOR_CHOICES = (
        (BILIRUBIN, "Bilirubin"),
        (BILIVERDIN, "Biliverdin"),
        (FLAVIN, "Flavin"),
        (PHYCOCYANOBILIN, "Phycocyanobilin"),
        (RIBITYL_LUMAZINE, "ribityl-lumazine"),
    )

    # Attributes
    # uuid        = models.UUIDField(default=uuid_lib.uuid4, editable=False, unique=True)  # for API
    uuid = models.CharField(
        max_length=5,
        default=prot_uuid,
        editable=False,
        unique=True,
        db_index=True,
        verbose_name="FPbase ID",
    )
    name = models.CharField(max_length=128, help_text="Name of the fluorescent protein", db_index=True)
    slug = models.SlugField(max_length=64, unique=True, help_text="URL slug for the protein")  # for generating urls
    base_name = models.CharField(max_length=128)  # easily searchable "family" name
    aliases = ArrayField(models.CharField(max_length=200), blank=True, null=True)
    chromophore = _NonNullChar(max_length=5, blank=True, default="")
    seq_validated = models.BooleanField(default=False, help_text="Sequence has been validated by a moderator")
    # seq must be nullable because of uniqueness contraints
    seq = SequenceField(
        unique=True,
        blank=True,
        null=True,
        verbose_name="Sequence",
        help_text="Amino acid sequence (IPG ID is preferred)",
    )
    seq_comment = models.CharField(
        max_length=512,
        blank=True,
        help_text="if necessary, comment on source of sequence",
    )

    pdb = ArrayField(
        models.CharField(max_length=4),
        blank=True,
        null=True,
        verbose_name="Protein DataBank IDs",
    )
    genbank = models.CharField(
        max_length=12,
        null=True,
        blank=True,
        unique=True,
        verbose_name="Genbank Accession",
        help_text="NCBI Genbank Accession",
    )
    uniprot = models.CharField(
        max_length=10,
        null=True,
        blank=True,
        unique=True,
        verbose_name="UniProtKB Accession",
        validators=[validate_uniprot],
    )
    ipg_id = models.CharField(
        max_length=12,
        null=True,
        blank=True,
        unique=True,
        verbose_name="IPG ID",
        help_text="Identical Protein Group ID at Pubmed",
    )  # identical protein group uid
    mw = models.FloatField(null=True, blank=True, help_text="Molecular Weight")  # molecular weight
    agg = models.CharField(
        max_length=2,
        choices=AGG_CHOICES,
        blank=True,
        verbose_name="Oligomerization",
        help_text="Oligomerization tendency",
    )
    oser = models.FloatField(null=True, blank=True, help_text="OSER score")  # molecular weight
    switch_type = models.CharField(
        max_length=2,
        choices=SWITCHING_CHOICES,
        blank=True,
        default="b",
        verbose_name="Switching Type",
        help_text="Photoswitching type (basic if none)",
    )
    blurb = models.TextField(max_length=512, blank=True, help_text="Brief descriptive blurb")
    cofactor = models.CharField(
        max_length=2,
        choices=COFACTOR_CHOICES,
        blank=True,
        help_text="Required for fluorescence",
    )

    # Relations
    parent_organism = models.ForeignKey(
        "Organism",
        related_name="proteins",
        verbose_name="Parental organism",
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
        help_text="Organism from which the protein was engineered",
    )
    primary_reference = models.ForeignKey(
        Reference,
        related_name="primary_proteins",
        verbose_name="Primary Reference",
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
        help_text="Preferably the publication that introduced the protein",
    )  # usually, the original paper that published the protein
    references = models.ManyToManyField(
        Reference, related_name="proteins", blank=True
    )  # all papers that reference the protein
    # FRET_partner = models.ManyToManyField('self', symmetrical=False, through='FRETpair', blank=True)
    default_state = models.ForeignKey(
        "State",
        related_name="default_for",
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )

    # __original_ipg_id = None

    # managers
    objects = ProteinManager()
    visible = QueryManager(~Q(status="hidden"))

    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     # store IPG_ID so that we know if it changes
    #     self.__original_ipg_id = self.ipg_id

    def mutations_from_root(self):
        try:
            root = self.lineage.get_root()
            if root.protein.seq and self.seq:
                return root.protein.seq.mutations_to(self.seq)
        except ObjectDoesNotExist:
            return None

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
        if self.status == "approved":
            return self
        try:
            return (
                Version.objects.get_for_object(self).filter(serialized_data__contains='"status": "approved"').first()
            )
        except Exception:
            return None

    @property
    def additional_references(self):
        return self.references.exclude(id=self.primary_reference_id).order_by("-year")

    @property
    def em_css(self):
        if self.states.count() > 1:
            from collections import OrderedDict

            stops = OrderedDict({st.emhex: "" for st in self.states.all()})
            bgs = []
            stepsize = int(100 / (len(stops) + 1))
            sub = 0
            for i, _hex in enumerate(stops):
                if _hex == "#000":
                    sub = 18
                bgs.append(f"{_hex} {(i + 1) * stepsize - sub}%")
            return f"linear-gradient(90deg, {', '.join(bgs)})"
        elif self.default_state:
            return self.default_state.emhex
        else:
            return "repeating-linear-gradient(-45deg,#333,#333 8px,#444 8px,#444 16px);"

    @property
    def em_svg(self):
        if self.states.count() <= 1:
            return self.default_state.emhex if self.default_state else "?"
        stops = [st.emhex for st in self.states.all()]
        stepsize = int(100 / (len(stops) + 1))
        svgdef = "linear:"
        for i, color in enumerate(stops):
            perc = (i + 1.0) * stepsize
            if color == "#000":
                perc *= 0.2
            svgdef += f'<stop offset="{perc}%" style="stop-color:{color};" />'
        return svgdef

    @property
    def color(self):
        try:
            return get_color_group(self.default_state.ex_max, self.default_state.em_max)[0]
        except Exception:
            return ""

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
        return self.seq.mutations_to(other, **kwargs) if self.seq and other else None

    def mutations_from(self, other, **kwargs):
        if isinstance(other, Protein):
            other = other.seq
        return other.seq.mutations_to(self.seq, **kwargs) if (self.seq and other) else None

    def has_spectra(self):
        return any(state.has_spectra() for state in self.states.all())

    def has_bleach_measurements(self):
        return self.states.filter(bleach_measurements__isnull=False).exists()

    def d3_spectra(self):
        spectra = []
        for state in self.states.all():
            spectra.extend(state.d3_dicts())
        return json.dumps(spectra)

    def spectra_img(self, fmt="svg", output=None, **kwargs):
        spectra = list(Spectrum.objects.filter(owner_state__protein=self).exclude(subtype="2p"))
        title = self.name if kwargs.pop("title", False) else None
        if kwargs.get("twitter", False):
            title = self.name
        info = ""
        if self.default_state:
            info += f"Ex/Em λ: {self.default_state.ex_max}/{self.default_state.em_max}"
            info += f"\nEC: {self.default_state.ext_coeff}  QY: {self.default_state.qy}"
        return spectra_fig(spectra, fmt, output, title=title, info=info, **kwargs)

    def set_default_state(self) -> bool:
        # FIXME: should allow control of default states in form
        # if only 1 state, make it the default state
        if not self.default_state or self.default_state.is_dark:
            if self.states.count() == 1 and not self.states.first().is_dark:
                self.default_state = self.states.first()
            # otherwise use farthest red non-dark state
            elif self.states.count() > 1:
                self.default_state = self.states.exclude(is_dark=True).order_by("-em_max").first()
            return True
        return False

    def clean(self):
        errors = {}
        # Don't allow basic switch_types to have more than one state.
        #        if self.switch_type == 'b' and self.states.count() > 1:
        #            errors.update({'switch_type': 'Basic (non photoconvertible) proteins '
        #                                          'cannot have more than one state.'})
        if self.pdb:
            self.pdb = list(set(self.pdb))
            for item in self.pdb:
                if Protein.objects.exclude(id=self.id).filter(pdb__contains=[item]).exists():
                    p = Protein.objects.filter(pdb__contains=[item]).first()
                    errors["pdb"] = f"PDB ID {item} is already in use by protein {p.name}"

        if errors:
            raise ValidationError(errors)

    def save(self, *args, **kwargs):
        # if the IPG ID has changed... refetch the sequence
        # if self.ipg_id != self.__original_ipg_id:
        #    s = fetch_ipg_sequence(uid=self.ipg_id)
        #    self.seq = s[1] if s else None

        self.slug = slugify(self.name)
        self.base_name = self._base_name
        super().save(*args, **kwargs)
        if self.set_default_state():
            super().save()
        # self.__original_ipg_id = self.ipg_id

    # Meta
    class Meta:
        ordering = ["name"]

    def history(self, ignoreKeys=()):
        from proteins.util.history import get_history

        return get_history(self, ignoreKeys)

    # ##################################
    # for algolia index

    def is_visible(self):
        return self.status != "hidden"

    def img_url(self):
        if self.has_spectra():
            return (
                "https://www.fpbase.org"
                + reverse("proteins:spectra-img", args=[self.slug])
                + ".png?xlabels=0&xlim=400,800"
            )
        else:
            return None

    def tags(self):
        tags = [self.get_switch_type_display(), self.get_agg_display(), self.color]
        return [i for i in tags if i]

    def date_published(self, norm=False):
        d = self.primary_reference.date if self.primary_reference else None
        if norm:
            return (d.year - 1992) / (datetime.datetime.now(pytz.utc).year - 1992) if d else 0
        return datetime.datetime.combine(d, datetime.datetime.min.time()) if d else None

    def n_faves(self, norm=False):
        nf = Favorite.objects.for_model(Protein).filter(target_object_id=self.id).count()
        if norm:
            from collections import Counter

            mx = Counter(Favorite.objects.for_model(Protein).values_list("target_object_id", flat=True)).most_common(1)
            mx = mx[0][1] if mx else 1
            return nf / mx
        return nf

    def n_cols(self):
        return ProteinCollection.objects.filter(proteins=self.id).count()

    def ga_views(self, period="month", norm=False):
        from proteins.extrest.ga import cached_ga_popular

        try:
            hits = cached_ga_popular()[period]
            return next(
                (
                    rating / max(list(zip(*hits))[2]) if norm else rating
                    for slug, _name, rating in hits
                    if slug == self.slug
                ),
                0,
            )
        except Exception:
            return 0

    def switchType(self):
        return self.get_switch_type_display()

    def _agg(self):
        return self.get_agg_display()

    def url(self):
        return self.get_absolute_url()

    def ex(self):
        if not self.states.exists():
            return None
        ex = [s.ex_max for s in self.states.all()]
        return ex[0] if len(ex) == 1 else ex

    def em(self):
        if not self.states.exists():
            return None
        em = [s.em_max for s in self.states.all()]
        return em[0] if len(em) == 1 else em

    def pka(self):
        if not self.states.exists():
            return None
        n = [s.pka for s in self.states.all()]
        return n[0] if len(n) == 1 else n

    def ec(self):
        if not self.states.exists():
            return None
        n = [s.ext_coeff for s in self.states.all()]
        return n[0] if len(n) == 1 else n

    def qy(self):
        if not self.states.exists():
            return None
        n = [s.qy for s in self.states.all()]
        return n[0] if len(n) == 1 else n

    def rank(self):
        # max rank is 1
        return (
            0.5 * self.date_published(norm=True) + 0.6 * self.ga_views(norm=True) + 1.0 * self.n_faves(norm=True)
        ) / 2.5

    def local_brightness(self):
        if self.states.exists():
            return max(s.local_brightness for s in self.states.all())

    def first_author(self):
        if self.primary_reference and self.primary_reference.first_author:
            return self.primary_reference.first_author.family
