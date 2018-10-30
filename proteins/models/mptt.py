from django.db import models
from mptt.models import MPTTModel, TreeForeignKey
from model_utils.models import TimeStampedModel
from references.models import Reference


class Lineage(MPTTModel, TimeStampedModel):
    protein = models.OneToOneField('Protein', on_delete=models.CASCADE, related_name='lineage')
    parent = TreeForeignKey('self', on_delete=models.CASCADE, null=True, blank=True, related_name='children')
    reference = models.ForeignKey(Reference, on_delete=models.CASCADE, null=True, blank=True, related_name='lineages')
    mutation = models.CharField(max_length=400, blank=True)

    class MPTTMeta:
        order_insertion_by = ['protein']

    def __repr__(self):
        return '<Lineage: {}>'.format(self)

    def __str__(self):
        return str(self.id)

    def display_mutation(self, maxwidth=30):
        if not self.mutation:
            return None
        if len(self.mutation) > maxwidth:
            muts = self.mutation.split('/')
            line = ''
            while len(line) < maxwidth - 4:
                if len(line):
                    line += '/'
                line += muts.pop(0)
            if len(line) > maxwidth - 4:
                line = line[:maxwidth] + '...'
            if len(muts):
                line += ' (+{})'.format(len(muts))
            return line
        return self.mutation
