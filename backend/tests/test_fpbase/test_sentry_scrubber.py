"""Password form fields are redacted from Sentry events."""

from __future__ import annotations

import pytest
from allauth.account import forms as account_forms
from django import forms
from django.conf import settings
from sentry_sdk.scrubber import EventScrubber

from fpbase.forms import CustomSignupForm

PASSWORD_FORMS = [
    CustomSignupForm,
    account_forms.LoginForm,
    account_forms.ChangePasswordForm,
    account_forms.SetPasswordForm,
    account_forms.ResetPasswordKeyForm,
]


@pytest.mark.django_db
@pytest.mark.parametrize("form_class", PASSWORD_FORMS, ids=lambda f: f.__name__)
def test_password_fields_are_scrubbed(form_class: type[forms.Form]):
    fields = form_class().fields
    passwords = [k for k, f in fields.items() if isinstance(f.widget, forms.PasswordInput)]
    assert passwords
    event = {"request": {"data": dict.fromkeys(fields, "value")}}
    EventScrubber(denylist=settings.SENTRY_DENYLIST, recursive=True).scrub_event(event)
    assert not [k for k in passwords if event["request"]["data"][k] == "value"]
