from django import forms
from django.contrib.auth.models import User

from references.models import Author, Reference


class UserMixin:
    def clean_user(self):
        if not self.cleaned_data["added_by"]:
            return User()
        return self.cleaned_data["added_by"]

    def clean_updated_by(self):
        if not self.cleaned_data["updated_by"]:
            return User()
        return self.cleaned_data["updated_by"]


class AuthorForm(UserMixin, forms.ModelForm):
    class Meta:
        model = Author
        fields = ["family", "given"]


class ReferenceForm(UserMixin, forms.ModelForm):
    refetch_info_on_save = forms.BooleanField(required=False)

    class Meta:
        model = Reference
        fields = ["pmid", "doi", "title", "refetch_info_on_save"]
