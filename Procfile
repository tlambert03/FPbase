release: python manage.py migrate --noinput
worker: celery worker --app=fpbase --concurrency 4 --without-gossip --without-mingle --without-heartbeat
web: gunicorn config.wsgi:application
