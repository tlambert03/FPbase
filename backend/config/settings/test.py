"""
Test settings for FPbase project.

- Used to run tests fast on the continuous integration server and locally
"""

import getpass
import os

from .base import *  # noqa

# DEBUG
# ------------------------------------------------------------------------------
# Turn debug off so tests run faster
DEBUG = False
TEMPLATES[0]["OPTIONS"]["debug"] = True

# ALLOWED_HOSTS
# ------------------------------------------------------------------------------
# Allow all hosts for testing (live_server uses random ports)
ALLOWED_HOSTS = ["*"]

# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This key only used for development and testing.
SECRET_KEY = env("DJANGO_SECRET_KEY", default="CHANGEME!!!")

# Use DATABASE_URL if set (e.g., in CI), otherwise use local test database
# Always use local database unless DATABASE_URL is explicitly valid (not placeholder values)
db_url = os.getenv("DATABASE_URL", "")
if not db_url or "user:password" in db_url:
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.postgresql",
            "USER": os.getenv("USER", getpass.getuser()),
            "PORT": "5432",
        }
    }

# Mail settings
# ------------------------------------------------------------------------------
EMAIL_HOST = "localhost"
EMAIL_PORT = 1025

# In-memory email backend stores messages in django.core.mail.outbox
# for unit testing purposes
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"

# CACHING
# ------------------------------------------------------------------------------
# Speed advantages of in-memory caching without having to run Memcached
CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.locmem.LocMemCache",
        "LOCATION": "",
    }
}

# CELERY
# ------------------------------------------------------------------------------
# Use eager mode for tests - tasks execute synchronously without Redis
CELERY_TASK_ALWAYS_EAGER = True
CELERY_TASK_EAGER_PROPAGATES = True
CELERY_TASK_STORE_EAGER_RESULT = True  # Store results even in eager mode for result.ready() checks
# Use in-memory broker and backend for tests
CELERY_BROKER_URL = "memory://"
CELERY_RESULT_BACKEND = "cache+memory://"

# TESTING
# ------------------------------------------------------------------------------
TEST_RUNNER = "django.test.runner.DiscoverRunner"
TESTING = True  # Flag for test-specific behavior

# RECAPTCHA
# ------------------------------------------------------------------------------
# Use Google's official test keys for reCAPTCHA testing
# See: https://developers.google.com/recaptcha/docs/faq#id-like-to-run-automated-tests-with-recaptcha.-what-should-i-do
RECAPTCHA_PUBLIC_KEY = "6LeIxAcTAAAAAJcZVRqyHh71UMIEGNQ_MXjiZKhI"
RECAPTCHA_PRIVATE_KEY = "6LeIxAcTAAAAAGG-vFI1TnRWxMZNFuojJ4WifJWe"
# Silence test key warnings
SILENCED_SYSTEM_CHECKS = ["django_recaptcha.recaptcha_test_key_error"]

# PASSWORD HASHING
# ------------------------------------------------------------------------------
# Use fast password hasher so tests run faster
PASSWORD_HASHERS = [
    "django.contrib.auth.hashers.MD5PasswordHasher",
]

# TEMPLATE LOADERS
# ------------------------------------------------------------------------------
# Keep templates in memory so tests run faster
TEMPLATES[0]["OPTIONS"]["loaders"] = [
    [
        "django.template.loaders.cached.Loader",
        [
            "django.template.loaders.filesystem.Loader",
            "django.template.loaders.app_directories.Loader",
        ],
    ],
]


# django-vite in test mode uses manifest (never dev server)
DJANGO_VITE = {
    "default": {
        "dev_mode": False,
        "manifest_path": str(ROOT_DIR.parent / "frontend" / "dist" / "manifest.json"),
    }
}
