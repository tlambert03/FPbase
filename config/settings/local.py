"""
Local settings for FPbase project.

- Run in Debug mode

- Use console backend for emails

- Add Django Debug Toolbar
- Add django-extensions as app
"""

from .base import *  # noqa


# DEBUG
# ------------------------------------------------------------------------------
DEBUG = env.bool('DJANGO_DEBUG', default=True)
TEMPLATES[0]['OPTIONS']['debug'] = DEBUG

CRISPY_FAIL_SILENTLY = not DEBUG

CSRF_COOKIE_HTTPONLY = True

# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This random key only used for development and testing, not on live site.
SECRET_KEY = env('DJANGO_SECRET_KEY', default='w)CU)uzJ<JMlkGTrfz?:)W>]EG!PFngIvQZq#9.r=sfHUmCPIe')

# Mail settings
# ------------------------------------------------------------------------------

EMAIL_PORT = 1025

EMAIL_HOST = 'localhost'
EMAIL_BACKEND = env('DJANGO_EMAIL_BACKEND',
                    default='django.core.mail.backends.console.EmailBackend')

if env('MAILGUN_API_KEY', default=False) and env('MAILGUN_DOMAIN', default=False):
    INSTALLED_APPS += ['anymail', ]
    ANYMAIL = {
        'MAILGUN_API_KEY': env('MAILGUN_API_KEY'),
        'MAILGUN_SENDER_DOMAIN': env('MAILGUN_DOMAIN')
    }
    EMAIL_BACKEND = 'anymail.backends.mailgun.EmailBackend'

ALLOWED_HOSTS = env.list('DJANGO_ALLOWED_HOSTS', default=['fpbase.org', 'localhost', 'testserver'])
# CORS
# -------

MIDDLEWARE = ['corsheaders.middleware.CorsMiddleware', ] + MIDDLEWARE
INSTALLED_APPS += ['corsheaders', ]
CORS_ORIGIN_WHITELIST = (
    'localhost:8000',
    'localhost:3000',
)
CORS_URLS_REGEX = r'^/test/.*$'


# CACHING
# ------------------------------------------------------------------------------
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': ''
    }
}

# django-debug-toolbar
# ------------------------------------------------------------------------------
MIDDLEWARE += ['debug_toolbar.middleware.DebugToolbarMiddleware', ]
INSTALLED_APPS += ['debug_toolbar', ]

INTERNAL_IPS = ['127.0.0.1', '10.0.2.2', ]

DEBUG_TOOLBAR_CONFIG = {
    'DISABLE_PANELS': [
        'debug_toolbar.panels.redirects.RedirectsPanel',
    ],
    'SHOW_TEMPLATE_CONTEXT': True,
}

# django-extensions
# ------------------------------------------------------------------------------
INSTALLED_APPS += ['django_extensions', ]

# TESTING
# ------------------------------------------------------------------------------
TEST_RUNNER = 'django.test.runner.DiscoverRunner'

# Your local stuff: Below this line define 3rd party library settings
# ------------------------------------------------------------------------------

#SITE_ID = None
