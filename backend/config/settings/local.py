"""
Local settings for FPbase project.

- Run in Debug mode

- Use console backend for emails

- Add Django Debug Toolbar
- Add django-extensions as app
"""

import structlog

from .base import *  # noqa

# STATIC FILES - Add backend static directory for development
# ------------------------------------------------------------------------------
# Include backend/fpbase/static so Django can serve microscope.js and other
# static files that don't go through webpack
STATICFILES_DIRS = [*STATICFILES_DIRS, str(ROOT_DIR / "fpbase" / "static")]

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = env.bool("DJANGO_DEBUG", default=True)
TEMPLATES[0]["OPTIONS"]["debug"] = DEBUG

CRISPY_FAIL_SILENTLY = not DEBUG

# DJANGO_VITE - Enable dev mode for local development
# ------------------------------------------------------------------------------
DJANGO_VITE["default"]["dev_mode"] = True

# CSRF_COOKIE_HTTPONLY = True

# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This random key only used for development and testing, not on live site.
SECRET_KEY = env("DJANGO_SECRET_KEY", default="w)CU)uzJ<JMlkGTrfz?:)W>]EG!PFngIvQZq#9.r=sfHUmCPIe")

# Mail settings
# ------------------------------------------------------------------------------

EMAIL_PORT = 1025

EMAIL_HOST = "localhost"
EMAIL_BACKEND = env("DJANGO_EMAIL_BACKEND", default="django.core.mail.backends.console.EmailBackend")

if env("MAILGUN_API_KEY", default=False) and env("MAILGUN_DOMAIN", default=False):
    INSTALLED_APPS += [
        "anymail",
    ]
    ANYMAIL = {
        "MAILGUN_API_KEY": env("MAILGUN_API_KEY"),
        "MAILGUN_SENDER_DOMAIN": env("MAILGUN_DOMAIN"),
    }
    EMAIL_BACKEND = "anymail.backends.mailgun.EmailBackend"

ALLOWED_HOSTS = env.list(
    "DJANGO_ALLOWED_HOSTS",
    default=["fpbase.org", "localhost", "testserver", "10.0.2.2", "127.0.0.1"],
)

# CACHING
# ------------------------------------------------------------------------------
CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.locmem.LocMemCache",
        "LOCATION": "",
    }
}


# django-debug-toolbar
# ------------------------------------------------------------------------------
MIDDLEWARE += [
    "debug_toolbar.middleware.DebugToolbarMiddleware",
]
INSTALLED_APPS += [
    "debug_toolbar",
]

INTERNAL_IPS = [
    "127.0.0.1",
    "10.0.2.2",
]

DEBUG_TOOLBAR_CONFIG = {
    "DISABLE_PANELS": [
        "debug_toolbar.panels.redirects.RedirectsPanel",
        "debug_toolbar.panels.profiling.ProfilingPanel",
    ],
    "SHOW_TEMPLATE_CONTEXT": True,
}

# django-extensions
# ------------------------------------------------------------------------------
INSTALLED_APPS += [
    "django_extensions",
]

# TESTING
# ------------------------------------------------------------------------------
TEST_RUNNER = "django.test.runner.DiscoverRunner"

# SITE_ID = None


SHELL_PLUS_POST_IMPORTS = [
    ("proteins.util.helpers", ("getprot", "getmut", "showalign")),
    ("proteins.util", ("maintain", "_local")),
    (
        "fpseq",
        ("FPSeq", "from_fpbase", "MutationSet", "get_mutations", "mutate_sequence"),
    ),
]

# Structlog Configuration for Local Development
structlog.configure(
    processors=[
        *STRUCTLOG_SHARED_PROCESSORS,
        structlog.stdlib.filter_by_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.StackInfoRenderer(),
        structlog.dev.set_exc_info,  # Dev-only: enhanced exception formatting
        structlog.processors.UnicodeDecoder(),
        structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
    ],
    logger_factory=structlog.stdlib.LoggerFactory(),
    cache_logger_on_first_use=True,
)

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "colored": {
            "()": structlog.stdlib.ProcessorFormatter,
            "processors": [
                structlog.stdlib.ProcessorFormatter.remove_processors_meta,
                structlog.dev.ConsoleRenderer(colors=True),
            ],
            "foreign_pre_chain": [
                *STRUCTLOG_SHARED_PROCESSORS,
                structlog.stdlib.ExtraAdder(),  # Merge extra={"key": "value"} from stdlib logging
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,  # Format tracebacks from stdlib logging
            ],
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "colored",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "django.server": {
            "level": "WARNING",  # Hide normal requests, use structlog instead
            "propagate": False,
        },
        "debug_toolbar": {
            "level": "WARNING",  # Hide debug toolbar noise
            "propagate": False,
        },
    },
}
