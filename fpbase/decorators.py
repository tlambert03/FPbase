from functools import wraps

import requests
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import REDIRECT_FIELD_NAME
from django.contrib.auth.decorators import login_required
from django.utils.decorators import available_attrs

default_message = "Please log in, in order to see the requested page."


def user_passes_test(test_func, message=default_message):
    """
    Decorator for views that checks that the user passes the given test,
    setting a message in case of no success. The test should be a callable
    that takes the user object and returns True if the user passes.
    """

    def decorator(view_func):
        @wraps(view_func, assigned=available_attrs(view_func))
        def _wrapped_view(request, *args, **kwargs):
            if not test_func(request.user):
                messages.error(request, message)
            return view_func(request, *args, **kwargs)

        return _wrapped_view

    return decorator


def login_required_message(function=None, message=default_message):
    """
    Decorator for views that checks that the user is logged in, redirecting
    to the log-in page if necessary.
    """
    actual_decorator = user_passes_test(lambda u: u.is_authenticated, message=message)
    if function:
        return actual_decorator(function)
    return actual_decorator


def login_required_message_and_redirect(
    function=None,
    redirect_field_name=REDIRECT_FIELD_NAME,
    login_url=None,
    message=default_message,
):

    if function:
        return login_required_message(
            login_required(function, redirect_field_name, login_url), message
        )

    return lambda deferred_function: login_required_message_and_redirect(
        deferred_function, redirect_field_name, login_url, message
    )


def check_recaptcha(view_func):
    @wraps(view_func)
    def wrap(request, *args, **kwargs):
        request.recaptcha_is_valid = None
        if request.method == "POST":
            recaptcha_response = request.POST.get("g-recaptcha-response")
            data = {
                "secret": settings.GOOGLE_RECAPTCHA_SECRET_KEY,
                "response": recaptcha_response,
            }
            r = requests.post(
                "https://www.google.com/recaptcha/api/siteverify", data=data
            )
            result = r.json()
            if result["success"]:
                request.recaptcha_is_valid = True
            else:
                request.recaptcha_is_valid = False
                messages.error(request, "Invalid reCAPTCHA. Please try again.")
        return view_func(request, *args, **kwargs)

    return wrap
