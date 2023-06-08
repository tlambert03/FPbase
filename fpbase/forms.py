from allauth.account.forms import SignupForm
from captcha.fields import ReCaptchaField
from captcha.widgets import ReCaptchaV3
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.conf import settings
from django.core.mail import EmailMessage
from django.utils.html import escape
from django.utils.safestring import mark_safe


class ContactForm(forms.Form):
    name = forms.CharField()
    email = forms.EmailField()
    message = forms.CharField(widget=forms.Textarea)
    captcha = ReCaptchaField(label="", widget=ReCaptchaV3(attrs={"required_score": 0.85}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit("submit", "Submit", css_class="btn btn-primary"))
        self.helper.form_method = "POST"

    def friendly_email(self):
        fullname = self.cleaned_data["name"]
        address = self.cleaned_data["email"]
        return mark_safe("%s <%s>") % (escape(fullname), escape(address))

    def send_email(self):
        EmailMessage(
            f'FPbase contact from {self.cleaned_data["name"]}',
            self.cleaned_data["message"],
            to=[a[1] for a in settings.ADMINS],
            headers={"Reply-To": self.friendly_email()},
        ).send()


class CustomSignupForm(SignupForm):
    captcha = ReCaptchaField(label="", widget=ReCaptchaV3(attrs={"required_score": 0.85}))
    field_order = ["username", "email", "password1", "password2", "captcha"]
