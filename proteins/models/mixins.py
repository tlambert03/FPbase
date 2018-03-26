from django.db import models
from django.contrib.auth import get_user_model
User = get_user_model()


class Authorable(models.Model):
    created_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_author')
    updated_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_modifier')

    class Meta:
        abstract = True
