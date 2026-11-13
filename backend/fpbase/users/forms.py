from __future__ import annotations

from allauth.account.adapter import get_adapter
from django import forms

from fpbase.users.models import User


class UserUpdateForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ["username", "name"]
        help_texts = {
            "username": "Used to log in. Links that contain your old username "
            "(e.g. to your public microscopes) will stop working."
        }

    def clean_username(self) -> str:
        username = self.cleaned_data["username"]
        # skip the "already taken" check if only the case changed (it's still us)
        shallow = username.lower() == self.instance.username.lower()
        return get_adapter().clean_username(username, shallow=shallow)


class UserDeleteForm(forms.Form):
    confirm = forms.CharField(label="Type your username to confirm")

    def __init__(self, *args, user: User, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.user = user

    def clean_confirm(self) -> str:
        confirm = self.cleaned_data["confirm"]
        if confirm != self.user.username:
            raise forms.ValidationError("That doesn't match your username.")
        return confirm
