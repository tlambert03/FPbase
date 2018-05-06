release: python manage.py migrate --noinput
web: newrelic-admin run-program gunicorn config.wsgi:application

