from allauth.account.models import EmailAddress
from avatar.templatetags.avatar_tags import avatar_url
from django import forms
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin as AuthUserAdmin
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.db.models import Count, Exists, OuterRef
from django.urls import reverse
from django.utils.safestring import mark_safe

from .models import User


class MyUserChangeForm(UserChangeForm):
    class Meta(UserChangeForm.Meta):
        model = User


class MyUserCreationForm(UserCreationForm):
    error_message = UserCreationForm.error_messages.update(
        {"duplicate_username": "This username has already been taken."}
    )

    class Meta(UserCreationForm.Meta):
        model = User

    def clean_username(self):
        username = self.cleaned_data["username"]
        try:
            User.objects.get(username=username)
        except User.DoesNotExist:
            return username
        raise forms.ValidationError(self.error_messages["duplicate_username"])


@admin.register(User)
class MyUserAdmin(AuthUserAdmin):
    form = MyUserChangeForm
    add_form = MyUserCreationForm
    fieldsets = (
        ("User Profile", {"fields": ("avatar", "name", "verified", "microscopes", "collections")}),
        *AuthUserAdmin.fieldsets,
    )
    list_display = (
        "username",
        "email",
        "_date_joined",
        "_last_login",
        "verified",
        "_logins",
        "social",
        "cols",
        "scopes",
        "faves",
    )
    search_fields = ["name", "username", "email"]
    readonly_fields = ("avatar", "verified", "social", "microscopes", "collections")
    ordering = ("-date_joined",)

    @admin.display(ordering="date_joined")
    def _date_joined(self, obj):
        return obj.date_joined.strftime("%Y/%m/%d")

    @admin.display(ordering="last_login")
    def _last_login(self, obj):
        if obj.last_login:
            return obj.last_login.strftime("%Y/%m/%d")
        return ""

    def avatar(self, obj):
        url = f'<img src="{avatar_url(obj)}" />'
        return mark_safe(url)

    @admin.display(boolean=True)
    def verified(self, obj):
        return obj.verified

    @admin.display(description="microscopes")
    def microscopes(self, obj):
        def _makelink(m):
            url = reverse("admin:proteins_microscope_change", args=(m.pk,))
            return f'<a href="{url}">{m.name}</a>'

        links = [_makelink(m) for m in obj.microscopes.all()]
        return mark_safe(", ".join(links))

    @admin.display(description="collections")
    def collections(self, obj):
        def _makelink(m):
            url = reverse("proteins:collection-detail", args=(m.pk,))
            return f'<a href="{url}">{m.name}</a>'

        links = [_makelink(m) for m in obj.proteincollections.all()]
        return mark_safe(", ".join(links))

    def social(self, obj):
        return ", ".join([q.provider.title() for q in obj.socialaccount_set.all()])

    @admin.display(ordering="_microscopes")
    def scopes(self, obj):
        return obj._microscopes or ""

    @admin.display(ordering="_favorites")
    def faves(self, obj):
        return obj._favorites or ""

    @admin.display(ordering="_logins")
    def _logins(self, obj):
        return obj._logins or ""

    @admin.display(ordering="_collections")
    def cols(self, obj):
        return obj._collections or ""

    def get_queryset(self, request):
        return (
            super()
            .get_queryset(request)
            .prefetch_related("socialaccount_set", "proteincollections", "emailaddress_set")
            .annotate(verified=Exists(EmailAddress.objects.filter(user_id=OuterRef("id"), verified=True)))
            .annotate(
                _collections=Count("proteincollections"),
                _microscopes=Count("microscopes"),
                _favorites=Count("favorites"),
                _logins=Count("logins"),
            )
        )
