from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING

from django.contrib import messages
from django.contrib.auth import REDIRECT_FIELD_NAME
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse

from fpbase.etag_utils import generate_version_etag, parse_etag_header

if TYPE_CHECKING:
    from collections.abc import Callable

    from django.db.models import Model
    from django.http import HttpRequest

default_message = "Please log in, in order to see the requested page."


def user_passes_test(test_func, message=default_message):
    """
    Decorator for views that checks that the user passes the given test,
    setting a message in case of no success. The test should be a callable
    that takes the user object and returns True if the user passes.
    """

    def decorator(view_func):
        @wraps(view_func)
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
        return login_required_message(login_required(function, redirect_field_name, login_url), message)

    return lambda deferred_function: login_required_message_and_redirect(
        deferred_function, redirect_field_name, login_url, message
    )


def etag_cached(*models: type[Model]) -> Callable:
    """Add ETag support to function-based views."""

    def decorator(view_func: Callable) -> Callable:
        @wraps(view_func)
        def wrapper(request: HttpRequest, *args, **kwargs) -> HttpResponse:
            current_etag = generate_version_etag(*models)

            if request.method in ("GET", "HEAD"):
                if_none_match = request.headers.get("if-none-match")
                if if_none_match:
                    client_etags = parse_etag_header(if_none_match)
                    if "*" in client_etags or current_etag in client_etags:
                        response = HttpResponse(status=304)
                        response["ETag"] = current_etag
                        return response

            response = view_func(request, *args, **kwargs)
            if response.status_code == 200:
                response["ETag"] = current_etag
            return response

        return wrapper

    return decorator
