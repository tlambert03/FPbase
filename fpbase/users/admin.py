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
    list_display = ('username', 'email', '_date_joined', '_last_login', '_logins', '_collections', 'social', )
    search_fields = ['name']
    readonly_fields = ('email_verified', 'social')

    def _date_joined(self, obj):
        return obj.date_joined.strftime("%Y/%m/%d")
    _date_joined.admin_order_field = 'date_joined'

    def _last_login(self, obj):
        if obj.last_login:
            return obj.last_login.strftime("%Y/%m/%d")
        return ''
    _last_login.admin_order_field = 'last_login'

    def email_verified(self, obj):
        return any([e.verified for e in obj.emailaddress_set.all()])
    email_verified.boolean = True

    def social(self, obj):
        return ", ".join([q.provider.title() for q in obj.socialaccount_set.all()])

    def _collections(self, obj):
        return obj._collections or ''
    _collections.admin_order_field = '_collections'

    def _logins(self, obj):
        return obj._logins or ''
    _logins.admin_order_field = '_logins'

    def get_queryset(self, request):
            return super(MyUserAdmin, self).get_queryset(request) \
                .prefetch_related('socialaccount_set',
                                  'proteincollections',
                                  'emailaddress_set') \
                .annotate(_collections=Count('proteincollections'),
                          _logins=Count('logins'))

