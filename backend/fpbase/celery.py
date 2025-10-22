import os

import scout_apm.celery
from celery import Celery
from django.conf import settings

if "DJANGO_SETTINGS_MODULE" not in os.environ:
    os.environ["DJANGO_SETTINGS_MODULE"] = "config.settings.local"

app = Celery("fpbase", namespace="CELERY")
app.config_from_object("django.conf:settings", namespace="CELERY")

scout_apm.celery.install(app)

# Load task modules from all registered Django app configs.
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)
