release: python manage.py migrate --noinput
worker: celery worker --app=fpbase --concurrency 4
web: gunicorn config.wsgi:application

