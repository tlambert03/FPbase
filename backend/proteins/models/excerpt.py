from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import models
from django.db.models import Q
from model_utils import Choices
from model_utils.managers import QueryManager
from model_utils.models import StatusModel, TimeStampedModel

from proteins.models.mixins import Authorable
from references.models import Reference

if TYPE_CHECKING:
    from proteins.models import Protein


class Excerpt(Authorable, TimeStampedModel, StatusModel):
    STATUS = Choices("approved", "flagged", "rejected")

    content = models.TextField(max_length=1200, help_text="Brief excerpt describing this protein")
    proteins: models.ManyToManyField[Protein, Excerpt] = models.ManyToManyField(
        "Protein", blank=True, related_name="excerpts"
    )
    reference_id: int | None
    reference: models.ForeignKey[Reference | None] = models.ForeignKey(
        Reference,
        related_name="excerpts",
        null=True,
        on_delete=models.SET_NULL,
        help_text="Source of this excerpt",
    )

    objects = models.Manager()
    visible = QueryManager(~Q(status="rejected"))

    class Meta:
        ordering = ["reference__year", "created"]

    def __str__(self):
        ref = self.reference.citation if self.reference else ""
        return f"{ref}: {self.content[:30]}..."

    def get_absolute_url(self):
        return self.reference.get_absolute_url()
