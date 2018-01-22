from django.contrib.auth.models import User
from django.forms import ModelForm
from references.models import Reference, Author


class UserMixin(object):
    def clean_user(self):
        if not self.cleaned_data['added_by']:
            return User()
        return self.cleaned_data['added_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']


class AuthorForm(UserMixin, ModelForm):
    class Meta:
        model = Author
        fields = ['family', 'given']


class ReferenceForm(UserMixin, ModelForm):
    class Meta:
        model = Reference
        fields = ['pmid', 'doi', ]


