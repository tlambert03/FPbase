from __future__ import annotations

from typing import TYPE_CHECKING

from django.contrib.auth.models import AbstractUser
from django.contrib.auth.signals import user_logged_in
from django.db import models
from django.urls import reverse
from django.utils.translation import gettext_lazy as _

if TYPE_CHECKING:
    from favit.models import Favorite
    from proteins.models import Microscope, OpticalConfig, ProteinCollection
    from references.models import Reference


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
    name = models.CharField(_("Name of User"), blank=True, max_length=255)

    if TYPE_CHECKING:
        logins: models.QuerySet[UserLogin]
        favorites: models.QuerySet[Favorite]
        reference_author: models.QuerySet[Reference]
        reference_modifier: models.QuerySet[Reference]
        microscopes: models.QuerySet[Microscope]
        opticalconfigs: models.QuerySet[OpticalConfig]
        proteincollections: models.QuerySet[ProteinCollection]

    def __str__(self):
        return self.username

    def get_absolute_url(self):
        return reverse("users:detail", kwargs={"username": self.username})


class UserLogin(models.Model):
    """Represent users' logins, one per record"""

    user_id: int
    user: models.ForeignKey[User] = models.ForeignKey(User, on_delete=models.CASCADE, related_name="logins")
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self) -> str:
        return f"{self.user} logged in at {self.timestamp}"


def update_user_login(sender, user, **kwargs):
    user.logins.create()
    user.save()


user_logged_in.connect(update_user_login)
