"""
ASGI entrypoint. Configures Django and then runs the application
defined in the ASGI_APPLICATION setting.
"""

import os
import django
from channels.routing import get_default_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.production")
django.setup()
application = get_default_application()
# if os.environ.get('DJANGO_SETTINGS_MODULE') == 'config.settings.production':
#     from raven.contrib.django.raven_compat.middleware.wsgi import Sentry
#     application = Sentry(application)
