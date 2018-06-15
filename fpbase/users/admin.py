from django import forms
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin as AuthUserAdmin
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.db.models import Count
from .models import User


class MyUserChangeForm(UserChangeForm):
    class Meta(UserChangeForm.Meta):
        model = User


class MyUserCreationForm(UserCreationForm):

    error_message = UserCreationForm.error_messages.update({
        'duplicate_username': 'This username has already been taken.'
    })

    class Meta(UserCreationForm.Meta):
        model = User

    def clean_username(self):
        username = self.cleaned_data["username"]
        try:
            User.objects.get(username=username)
        except User.DoesNotExist:
            return username
        raise forms.ValidationError(self.error_messages['duplicate_username'])


@admin.register(User)
class MyUserAdmin(AuthUserAdmin):
    form = MyUserChangeForm
    add_form = MyUserCreationForm
    fieldsets = (
        ('User Profile', {'fields': ('name', 'email_verified')}),
    ) + AuthUserAdmin.fieldsets
    list_display = ('username', 'email', 'join_date', 'email_verified', 'social', '_collections')
    search_fields = ['name']
    readonly_fields = ('email_verified', 'social')

    def join_date(self, obj):
        return obj.date_joined.strftime("%Y/%m/%d")
    join_date.admin_order_field = 'date_joined'

    def email_verified(self, obj):
        return any([e.verified for e in obj.emailaddress_set.all()])
    email_verified.boolean = True

    def social(self, obj):
        return ", ".join([q.provider.title() for q in obj.socialaccount_set.all()])

    def _collections(self, obj):
        return obj._collections or ''
    _collections.admin_order_field = '_collections'

    def get_queryset(self, request):
            return super(MyUserAdmin, self).get_queryset(request) \
                .prefetch_related('socialaccount_set',
                                  'collections',
                                  'emailaddress_set') \
                .annotate(_collections=Count('collections'))

