from __future__ import annotations

from typing import TYPE_CHECKING

from allauth.account.adapter import get_adapter
from allauth.account.decorators import reauthentication_required
from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.messages.views import SuccessMessageMixin
from django.http import Http404
from django.shortcuts import redirect, render
from django.urls import reverse
from django.views.generic import DetailView, ListView, RedirectView, UpdateView

from fpbase.users.forms import UserDeleteForm, UserUpdateForm
from fpbase.users.models import User

if TYPE_CHECKING:
    from django.http import HttpRequest, HttpResponse


class UserDetailView(LoginRequiredMixin, DetailView):
    model = User
    # These next two lines tell the view to index lookups by username
    slug_field = "username"
    slug_url_kwarg = "username"

    def get_object(self, queryset=None):
        # a profile is visible to its owner (and to staff, from the user list)
        obj = super().get_object(queryset)
        if obj != self.request.user and not self.request.user.is_staff:
            raise Http404
        return obj


class UserRedirectView(LoginRequiredMixin, RedirectView):
    permanent = False

    def get_redirect_url(self):
        return reverse("users:detail", kwargs={"username": self.request.user.username})


class UserUpdateView(LoginRequiredMixin, SuccessMessageMixin, UpdateView):
    form_class = UserUpdateForm
    success_message = "Your profile was updated."

    def get_success_url(self):
        return reverse("users:detail", kwargs={"username": self.object.username})

    def get_object(self, queryset=None):
        # a fresh copy, so an invalid form doesn't mutate request.user
        return User.objects.get(pk=self.request.user.pk)

    def form_valid(self, form):
        response = super().form_valid(form)
        if "username" in form.changed_data:
            get_adapter(self.request).send_notification_mail(
                "account/email/username_changed",
                self.object,
                {"old_username": form.initial["username"], "new_username": self.object.username},
            )
        return response


class UserListView(LoginRequiredMixin, ListView):
    model = User
    # These next two lines tell the view to index lookups by username
    slug_field = "username"
    slug_url_kwarg = "username"


@login_required
@reauthentication_required
def user_delete(request: HttpRequest) -> HttpResponse:
    user = request.user
    if user.is_staff:
        messages.error(request, "Staff accounts must be deleted by an admin.")
        return redirect("users:detail", username=user.username)

    form = UserDeleteForm(request.POST if request.method == "POST" else None, user=user)
    if form.is_valid():
        # sent first, while the user's email addresses still exist (ATOMIC_REQUESTS
        # rolls back the delete if sending fails)
        get_adapter(request).send_notification_mail(
            "account/email/account_deleted", user, {"username": user.username}
        )
        # collections go with the user; everything else they contributed
        # (microscopes, proteins, spectra, references...) stays, unattributed.
        user.proteincollections.all().delete()
        user.delete()
        logout(request)
        messages.success(request, "Your account has been deleted.")
        return redirect("/")
    return render(request, "users/user_confirm_delete.html", {"form": form})
