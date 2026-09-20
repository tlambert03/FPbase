# Prepended to every remote script by remote.py.  Runs on a Heroku one-off dyno via
# `python -` from /app, so mirror what backend/manage.py does to sys.path.
# (backend/fpbase must be *appended*: it contains celery.py, which would shadow celery)
import json
import sys

sys.path.insert(0, "backend")
sys.path.append("backend/fpbase")

import django

django.setup()

from django.db import connection  # noqa: E402


def read_only():
    """Make Postgres refuse every write for the rest of this session."""
    connection.cursor().execute("SET default_transaction_read_only = on")


def emit(obj):
    """Print the result between markers so remote.py can ignore dyno log noise."""
    print("<<<JSON")
    print(json.dumps(obj, default=str))
    print("JSON>>>")
