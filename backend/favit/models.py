from __future__ import annotations

from typing import TYPE_CHECKING

from django.conf import settings
from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType
from django.db import models
from django.utils.translation import gettext_lazy as _

from .managers import FavoriteManager

if TYPE_CHECKING:
    from fpbase.users.models import User


class Favorite(models.Model):
    user_id: int
    user: models.ForeignKey[User] = models.ForeignKey(
        getattr(settings, "AUTH_USER_MODEL", "auth.User"),
        related_name="favorites",
        on_delete=models.CASCADE,
    )
    target_content_type_id: int
    target_content_type: models.ForeignKey[ContentType] = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    target_object_id = models.PositiveIntegerField()
    target = GenericForeignKey("target_content_type", "target_object_id")
    timestamp = models.DateTimeField(auto_now_add=True, db_index=True)

    objects: FavoriteManager = FavoriteManager()

    class Meta:
        ordering = ["-timestamp"]
        get_latest_by = "timestamp"
        unique_together = ("user", "target_content_type", "target_object_id")
        verbose_name = _("favorite")
        verbose_name_plural = _("favorites")

    def __str__(self) -> str:
        return f"{self.user} favorited {self.target}"
