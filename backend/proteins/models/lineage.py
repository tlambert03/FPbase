from django.core.exceptions import ValidationError
from django.db import models
from django.utils.translation import gettext_lazy as _
from model_utils.models import TimeStampedModel
from mptt.models import MPTTModel, TreeForeignKey

from fpseq.mutations import MutationSet
from references.models import Reference

from ..models.mixins import Authorable
from ..util.maintain import validate_node
from .protein import Protein


def parse_mutation(mut_string):
    try:
        return MutationSet(mut_string)
    except ValueError as e:
        raise ValidationError(_(f"Invalid input for MutationSet: {e}")) from e


class MutationSetField(models.CharField):
    description = _("String representing a protein mutation (up to %(max_length)s)")

    def from_db_value(self, value, expression, connection):
        return value if value is None else parse_mutation(value)

    def to_python(self, value):
        if isinstance(value, MutationSet):
            return value
        return value if value is None else parse_mutation(value)

    def get_prep_value(self, value):
        return str(value) if value else ""


class Lineage(MPTTModel, TimeStampedModel, Authorable):
    protein = models.OneToOneField("Protein", on_delete=models.CASCADE, related_name="lineage")
    parent = TreeForeignKey("self", on_delete=models.CASCADE, null=True, blank=True, related_name="children")
    reference = models.ForeignKey(
        Reference,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        related_name="lineages",
    )
    mutation = MutationSetField(max_length=400, blank=True)
    rootmut = models.CharField(max_length=400, blank=True)
    root_node = models.ForeignKey(
        "self",
        null=True,
        on_delete=models.CASCADE,
        related_name="descendants",
        verbose_name="Root Node",
    )

    class MPTTMeta:
        order_insertion_by = ["protein"]

    def save(self, *args, **kwargs):
        if not self.pk:
            kwargs["force_insert"] = False
            super().save(*args, **kwargs)
        if self.parent:
            try:
                self.root_node = self.get_root()
            except models.ObjectDoesNotExist:
                self.root_node = None
        else:
            self.root_node = None

        if self.pk and self.parent and self.parent.protein.seq:
            self.rootmut = self.mut_from_root()
        else:
            self.rootmut = ""
        super().save(*args, **kwargs)

    def mut_from_root(self, root=None):
        if root:
            if not isinstance(root, Protein):
                raise ValueError("root argument must be a protein instance")
        else:
            root = self.root_node.protein if self.root_node else self.get_root().protein
        # if the parent is the root, just return the mutation
        if root == self.parent.protein:
            return str(self.mutation)
        if not isinstance(self.mutation, MutationSet):
            self.mutation = parse_mutation(self.mutation)
        return self.mutation.relative_to_root(self.parent.protein.seq, root.seq)

    def __repr__(self):
        return f"<Lineage: {self}>"

    def __str__(self):
        return str(self.protein)

    def clean(self):
        errors = validate_node(self)
        if errors:
            raise ValidationError({"mutation": ValidationError(errors[0])})

    def derive_mutation(self, root=None):
        ms = None
        if self.parent and self.parent.protein.seq and self.protein.seq:
            if root:
                if not isinstance(root, Protein):
                    root = self.get_root().protein
                if root.seq:
                    return self.parent.protein.seq.mutations_to(self.protein.seq, reference=root.seq)
            else:
                ms = self.parent.protein.seq.mutations_to(self.protein.seq)
        return ms
