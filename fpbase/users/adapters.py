import dns.exception
import dns.resolver
from allauth.account.adapter import DefaultAccountAdapter
from allauth.socialaccount.adapter import DefaultSocialAccountAdapter
from django import forms
from django.conf import settings
from django.utils.translation import gettext_lazy as _


class AccountAdapter(DefaultAccountAdapter):
    def is_open_for_signup(self, request):
        return getattr(settings, "ACCOUNT_ALLOW_REGISTRATION", True)

    def clean_email(self, email):
        domain = email.split("@")[-1]
        try:
            dns.resolver.query(domain, "MX")
        except dns.exception.DNSException as e:
            raise forms.ValidationError(_("The domain %s could not be found.") % domain) from e
        return email


class SocialAccountAdapter(DefaultSocialAccountAdapter):
    def is_open_for_signup(self, request, sociallogin):
        return getattr(settings, "ACCOUNT_ALLOW_REGISTRATION", True)
