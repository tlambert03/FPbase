from django.db import models
from django.contrib.auth import get_user_model
User = get_user_model()


class Authorable(models.Model):
    created_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_author', on_delete=models.SET_NULL)
    updated_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_modifier', on_delete=models.SET_NULL)

    class Meta:
        abstract = True
