from django import forms
from django.core.mail import mail_admins


class ContactForm(forms.Form):
    name = forms.CharField()
    email = forms.EmailField()
    message = forms.CharField(widget=forms.Textarea)

    def send_email(self):
        print('hi')
        # send email using the self.cleaned_data dictionary
        mail_admins(
            'FPbase ContactForm Submission',
            'From: {}\n\n'.format(self.cleaned_data['name'])
            + 'Email: {}\n\n'.format(self.cleaned_data['email'])
            + self.cleaned_data['message'],
        )
