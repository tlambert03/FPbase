import os

from celery import Celery
from django.conf import settings

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.local")

app = Celery("fpbase", namespace="CELERY")
app.config_from_object("django.conf:settings", namespace="CELERY")

# Configure SSL settings for Redis if using secure connection
if settings.CELERY_RESULT_BACKEND.startswith("rediss://"):
    app.conf.update(
        redis_backend_use_ssl=settings.CELERY_REDIS_BACKEND_USE_SSL, broker_use_ssl=settings.CELERY_BROKER_USE_SSL
    )

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()
