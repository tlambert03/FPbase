"""
Local production-like settings for testing.

This file is used by `just prod-local` to test production configuration locally
without hitting production resources.

Differences from production:
- Uses local database, Redis, and dummy external services
- Redis SSL settings are disabled (local Redis doesn't use SSL)
- Celery SSL settings are disabled
- All external service credentials are dummy values
- Everything else matches production exactly
"""

import os

# Set dummy environment variables for local testing BEFORE importing production
os.environ.setdefault("DJANGO_SECRET_KEY", "local-production-test-key-change-in-real-production")
os.environ.setdefault("DATABASE_URL", "postgres:///fpbase")
os.environ.setdefault("DJANGO_ALLOWED_HOSTS", "localhost,127.0.0.1")
os.environ.setdefault("DJANGO_SECURE_SSL_REDIRECT", "False")
os.environ.setdefault("DJANGO_AWS_ACCESS_KEY_ID", "dummy-aws-key")
os.environ.setdefault("DJANGO_AWS_SECRET_ACCESS_KEY", "dummy-aws-secret")
os.environ.setdefault("MAILGUN_API_KEY", "key-888888888888888888")
os.environ.setdefault("MAILGUN_DOMAIN", "mg.fpbase.org")
os.environ.setdefault("SENTRY_DSN", "https://dummy@sentry.io/0")
os.environ.setdefault("SENTRY_PROJECT", "local-test")
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/0")
os.environ.setdefault("ALGOLIA_API_KEY", "")  # Empty to disable Algolia
os.environ.setdefault("SCOUT_MONITOR", "False")

from .production import *  # noqa

# Override Redis cache settings to disable SSL (local Redis doesn't use SSL)
CACHES = {
    "default": {
        "BACKEND": "django_redis.cache.RedisCache",
        "LOCATION": REDIS_LOCATION,
        "OPTIONS": {
            "CLIENT_CLASS": "django_redis.client.DefaultClient",
            "IGNORE_EXCEPTIONS": True,  # mimics memcache behavior.
            # No CONNECTION_POOL_KWARGS with ssl_cert_reqs for local testing
        },
    }
}

# Remove Celery SSL settings for local testing
g = globals()
g.pop("CELERY_BROKER_TRANSPORT_OPTIONS", None)
g.pop("CELERY_RESULT_BACKEND_TRANSPORT_OPTIONS", None)

# Disable secure cookie settings for localhost (HTTP, not HTTPS)
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False
SECURE_SSL_REDIRECT = False
