from allauth.account.forms import SignupForm
from captcha.fields import ReCaptchaField
from captcha.widgets import ReCaptchaV3
from django import forms
from django.conf import settings
from django.core.mail import EmailMessage
from django.utils.html import escape
from django.utils.safestring import mark_safe

# name = forms.CharField(widget=forms.TextInput(attrs={"class": "form-control"}))
# email = forms.EmailField(widget=forms.EmailInput(attrs={"class": "form-control"}))
# message = forms.CharField(widget=forms.Textarea(attrs={"class": "form-control"}))


class ContactForm(forms.Form):
    name = forms.CharField()
    email = forms.EmailField()
    message = forms.CharField(widget=forms.Textarea)
    captcha = ReCaptchaField(label="", widget=ReCaptchaV3(attrs={"required_score": 0.85}))

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
    captcha = ReCaptchaField(label="")
    field_order = ["username", "email", "password1", "password2", "captcha"]
