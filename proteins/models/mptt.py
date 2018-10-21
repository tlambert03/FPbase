from django.db import models
from mptt.models import MPTTModel, TreeForeignKey
from model_utils.models import TimeStampedModel


class Lineage(MPTTModel, TimeStampedModel):
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE, related_name='lineages')
    parent = TreeForeignKey('self', on_delete=models.CASCADE, null=True, blank=True, related_name='children')
    mutation = models.CharField(max_length=400, blank=True)

    class MPTTMeta:
        order_insertion_by = ['protein']

    def __repr__(self):
        return self.protein.name

    def __str__(self):
        return self.protein.name

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
