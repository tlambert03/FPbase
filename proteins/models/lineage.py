from django.db import models
from mptt.models import MPTTModel, TreeForeignKey
from model_utils.models import TimeStampedModel
from references.models import Reference
from django.utils.translation import gettext_lazy as _
from fpseq.mutations import MutationSet
from django.core.exceptions import ValidationError
from ..validators import validate_mutationset
from ..models.mixins import Authorable
from ..models import Protein


def parse_mutation(mut_string):
    validate_mutationset(mut_string)
    try:
        return MutationSet(mut_string)
    except ValueError as e:
        raise ValidationError(_("Invalid input for MutationSet: {}".format(e)))


class MutationSetField(models.CharField):
    description = _("String representing a protein mutation (up to %(max_length)s)")

    def from_db_value(self, value, expression, connection):
        if value is None:
            return value
        return parse_mutation(value)

    def to_python(self, value):
        if isinstance(value, MutationSet):
            return value
        if value is None:
            return value
        return parse_mutation(value)

    def get_prep_value(self, value):
        return str(value) if value else ''


class Lineage(MPTTModel, TimeStampedModel, Authorable):
    protein = models.OneToOneField('Protein', on_delete=models.CASCADE, related_name='lineage')
    parent = TreeForeignKey('self', on_delete=models.CASCADE, null=True, blank=True, related_name='children')
    reference = models.ForeignKey(Reference, on_delete=models.CASCADE, null=True, blank=True, related_name='lineages')
    mutation = MutationSetField(max_length=400, blank=True)
    rootmut = models.CharField(max_length=400, blank=True)

    class MPTTMeta:
        order_insertion_by = ['protein']

    def save(self, *args, **kwargs):
        if self.pk and self.parent and self.parent.seq and self.mutation:
            self.rootmut = self.mut_from_root()
        super().save(*args, **kwargs)

    def __repr__(self):
        return '<Lineage: {}>'.format(self)

    def __str__(self):
        return str(self.id)

    def mut_from_root(self, root=None):
        if root:
            if not isinstance(root, Protein):
                raise ValueError('root argument must be a protein instance')
        else:
            root = self.get_root().protein
        return self.mutation.relative_to_root(self.parent.protein.seq, root.seq)

    def derive_mutation(self, root=None):
        ms = None
        if self.parent and self.parent.protein.seq and self.protein.seq:
            if root:
                if not isinstance(root, Protein):
                    root = self.get_root().protein
                if root.seq:
                    ms = self.parent.protein.seq.mutations_to(self.protein.seq, reference=root.seq)
                    return ms
            else:
                ms = self.parent.protein.seq.mutations_to(self.protein.seq)
        return ms