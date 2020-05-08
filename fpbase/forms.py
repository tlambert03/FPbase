from allauth.account.forms import SignupForm
from captcha.fields import ReCaptchaField
from django import forms
from django.conf import settings
from django.core.mail import EmailMessage
from django.utils.html import escape
from django.utils.safestring import mark_safe


class ContactForm(forms.Form):
    name = forms.CharField()
    email = forms.EmailField()
    message = forms.CharField(widget=forms.Textarea)
    captcha = ReCaptchaField(label="")

    def friendly_email(self):
        fullname = self.cleaned_data["name"]
        address = self.cleaned_data["email"]
        return mark_safe(u"%s <%s>") % (escape(fullname), escape(address))

    def send_email(self):
        EmailMessage(
            "FPbase contact from %s" % self.cleaned_data["name"],
            self.cleaned_data["message"],
            # from_email=self.friendly_email(),
            to=[a[1] for a in settings.ADMINS],
            headers={"Reply-To": self.friendly_email()},
        ).send()


class CustomSignupForm(SignupForm):
    captcha = ReCaptchaField(label="")
    field_order = ["username", "email", "password1", "password2", "captcha"]
