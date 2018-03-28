from django.contrib.auth.models import AbstractUser
from django.urls import reverse
from django.db import models
from django.utils.translation import ugettext_lazy as _


class User(AbstractUser):

    # # AbstractUser Fields
    # username = models.CharField
    # first_name = models.CharField
    # last_name = models.CharField
    # email = models.EmailField
    # is_staff = models.BooleanField
    # is_active = models.BooleanField
    # date_joined = models.DateTimeField

    # First Name and Last Name do not cover name patterns
    # around the globe.
    name = models.CharField(_('Name of User'), blank=True, max_length=255)

    def __str__(self):
        return self.username

    def get_absolute_url(self):
        return reverse('users:detail', kwargs={'username': self.username})
