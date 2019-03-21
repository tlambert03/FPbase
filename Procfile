release: python manage.py migrate --noinput
worker: celery worker --app=fpbase -l info --concurrency 2
web: gunicorn config.wsgi:application

