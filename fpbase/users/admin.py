from django import forms
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin as AuthUserAdmin
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.db.models import Count, Exists, OuterRef
from .models import User
from avatar.templatetags.avatar_tags import avatar_url
from allauth.account.models import EmailAddress


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
        (
            "User Profile",
            {"fields": ("avatar", "name", "verified", "microscopes", "collections")},
        ),
    ) + AuthUserAdmin.fieldsets
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

    def _date_joined(self, obj):
        return obj.date_joined.strftime("%Y/%m/%d")

    _date_joined.admin_order_field = "date_joined"

    def _last_login(self, obj):
        if obj.last_login:
            return obj.last_login.strftime("%Y/%m/%d")
        return ""

    _last_login.admin_order_field = "last_login"

    def avatar(self, obj):
        url = '<img src="{}" />'.format(avatar_url(obj))
        return mark_safe(url)

    avatar.allow_tags = True

    def verified(self, obj):
        return obj.verified

    verified.boolean = True

    def microscopes(self, obj):
        def _makelink(m):
            url = reverse("admin:proteins_microscope_change", args=(m.pk,))
            return '<a href="{}">{}</a>'.format(url, m.name)

        links = [_makelink(m) for m in obj.microscopes.all()]
        return mark_safe(", ".join(links))

    microscopes.short_description = "microscopes"

    def collections(self, obj):
        def _makelink(m):
            url = reverse("proteins:collection-detail", args=(m.pk,))
            return '<a href="{}">{}</a>'.format(url, m.name)

        links = [_makelink(m) for m in obj.proteincollections.all()]
        return mark_safe(", ".join(links))

    collections.short_description = "collections"

    def social(self, obj):
        return ", ".join([q.provider.title() for q in obj.socialaccount_set.all()])

    def scopes(self, obj):
        return obj._microscopes or ""

    scopes.admin_order_field = "_microscopes"

    def faves(self, obj):
        return obj._favorites or ""

    faves.admin_order_field = "_favorites"

    def _logins(self, obj):
        return obj._logins or ""

    _logins.admin_order_field = "_logins"

    def cols(self, obj):
        return obj._collections or ""

    cols.admin_order_field = "_collections"

    def get_queryset(self, request):
        return (
            super(MyUserAdmin, self)
            .get_queryset(request)
            .prefetch_related(
                "socialaccount_set", "proteincollections", "emailaddress_set"
            )
            .annotate(
                verified=Exists(
                    EmailAddress.objects.filter(user_id=OuterRef("id"), verified=True)
                )
            )
            .annotate(
                _collections=Count("proteincollections"),
                _microscopes=Count("microscopes"),
                _favorites=Count("favorites"),
                _logins=Count("logins"),
            )
        )
