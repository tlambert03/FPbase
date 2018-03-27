# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import ValidationError
import re
from ..validators import validate_mutation

# WORK IN PROGRESS


class Mutations(models.Model):
    parent      = models.ForeignKey('Protein', related_name='proteins', verbose_name="Parent Protein")
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


# class Mutation(object):

#     def __init__(self, mutstring, delim=':'):
#         self.delim = delim
#         try:
#             splits = mutstring.split(delim)
#             parent = splits[0]
#             self.parent = Protein.object.get(id=parent)
#             mutations = splits[1]
#         except IndexError:
#             raise
#         except Exception:
#             raise

#         if mutations:
#             if isinstance(mutations, (list, set, tuple)):  # must be a list
#                 self.mutations = set(mutations)
#             elif isinstance(mutations, str):
#                 splits = mutations.split('/')
#                 for item in splits:
#                     validate_mutation(item)
#                 self.mutations = set(splits)
#             else:
#                 raise TypeError("mutations input must be list, set, tuple, or str")
#         else:
#             self.mutations = set()

#     def mutnum(self, val):
#         num = re.search(r'\d+', val)
#         if num:
#             return int(num.group())
#         else:
#             return None

#     def __str__(self):
#         ordered_mutations = list(self.mutations)
#         ordered_mutations.sort(key=self.mutnum)
#         return "{}{}{}".format(self.parent, self.delim, "/".join(ordered_mutations))


# class MutationsField(models.TextField):
#     description = "Stores a mutation object"

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#     def from_db_value(self, value, expression, connection, context):
#         if not value:
#             return None
#         return Mutation(value)

#     def to_python(self, value):
#         if isinstance(value, Mutation):
#             return value

#         if not value:
#             return None

#         try:
#             return Mutation(value)
#         except Exception:
#             raise ValidationError("Invalid input for a Mutation instance")

#     def get_prep_value(self, value):
#         if value is None:
#             return value
#         return str(value)
