from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING

from django.contrib import messages
from django.contrib.auth import REDIRECT_FIELD_NAME
from django.contrib.auth.decorators import login_required

from fpbase.etag_utils import check_etag_match, generate_version_etag

if TYPE_CHECKING:
    from collections.abc import Callable

    from django.db.models import Model
    from django.http import HttpRequest, HttpResponse


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


def etag_cached(
    *models: type[Model],
) -> Callable[[Callable[[HttpRequest], HttpResponse]], Callable[[HttpRequest], HttpResponse]]:
    """Add ETag support to function-based views."""

    def decorator(
        view_func: Callable[[HttpRequest], HttpResponse],
    ) -> Callable[[HttpRequest], HttpResponse]:
        @wraps(view_func)
        def wrapper(request: HttpRequest) -> HttpResponse:
            # Check if client's ETag matches - return 304 if so
            not_modified = check_etag_match(request, *models)
            if not_modified:
                return not_modified

            # Process the request normally
            response = view_func(request)

            # Add ETag header to successful responses
            if response.status_code == 200:
                response["ETag"] = generate_version_etag(*models)

            return response

        return wrapper

    return decorator
