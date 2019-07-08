# -*- coding: utf-8 -*-
from django.contrib import admin

from .models import Favorite


class FavoriteAdmin(admin.ModelAdmin):
    list_display = ["user", "target", "timestamp"]
    list_select_related = True
    search_fields = ("user__username",)
    raw_id_fields = ("user",)
    fields = ("user", "target_content_type", "target")
    readonly_fields = ("target",)

    def target(self, obj):
        cls = obj.target_content_type.model_class()
        obj = cls.objects.get(id=obj.target_object_id)
        return obj or ""

    target.admin_order_field = "target_object_id"


admin.site.register(Favorite, FavoriteAdmin)
