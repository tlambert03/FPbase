from django.db import models
from django.db.models import Q
from references.models import Reference
from model_utils.models import TimeStampedModel, StatusModel
from model_utils import Choices
from model_utils.managers import QueryManager
from .mixins import Authorable


class Excerpt(Authorable, TimeStampedModel, StatusModel):

    STATUS = Choices('approved', 'flagged', 'rejected')

    content = models.TextField(max_length=1024, help_text="Brief excerpt describing this protein")
    protein = models.ForeignKey('Protein', related_name='excerpts', on_delete=models.CASCADE)
    reference = models.ForeignKey(Reference, related_name='excerpt',
                                  null=True, on_delete=models.SET_NULL,
                                  help_text="Source of this excerpt")

    objects = models.Manager()
    visible = QueryManager(~Q(status='rejected'))

    class Meta:
        ordering = ['reference__year', 'created']

    def __str__(self):
        ref = self.reference.citation if self.reference else ''
        return "{}: {}...".format(ref, self.content[:30])

    def get_absolute_url(self):
        return self.protein.get_absolute_url()
