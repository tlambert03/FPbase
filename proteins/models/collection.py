from django.db import models
from django.contrib.auth import get_user_model
from django.urls import reverse
from model_utils.models import TimeStampedModel

User = get_user_model()


class ProteinCollection(TimeStampedModel):
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=512, blank=True)
    proteins = models.ManyToManyField('Protein', related_name='collection_memberships')
    owner = models.ForeignKey(User, blank=True, null=True, related_name='collections', verbose_name='Protein Collection', on_delete=models.SET_NULL,)
    private = models.BooleanField(default=False, verbose_name="Private Collection", help_text="Private collections can not be seen by or shared with other users")

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse("proteins:collection-detail", args=[self.id])

    class Meta:
        unique_together = (('owner', 'name'),)
