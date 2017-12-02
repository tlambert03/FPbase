from django.contrib.auth.models import User
from django.forms import ModelForm
from references.models import Reference, Author


class AuthorForm(ModelForm):
    class Meta:
        model = Author
        fields = ['last_name', 'initials']

    def clean_user(self):
        if not self.cleaned_data['added_by']:
            return User()
        return self.cleaned_data['added_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']


class ReferenceForm(ModelForm):
    class Meta:
        model = Reference
        fields = ['pmid', 'doi', ]

    def clean_user(self):
        if not self.cleaned_data['added_by']:
            return User()
        return self.cleaned_data['added_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']
