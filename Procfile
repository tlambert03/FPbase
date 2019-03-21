release: python manage.py migrate --noinput
worker: celery worker --app=fpbase
web: gunicorn config.wsgi:application

