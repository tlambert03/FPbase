# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import ValidationError
import re
from ..validators import validate_mutation


class Mutation(models.Model):
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
