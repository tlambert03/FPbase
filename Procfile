release: python manage.py migrate --noinput
worker: celery --app fpbase worker --concurrency 4 --without-gossip --without-mingle --without-heartbeat
web: gunicorn config.wsgi:application
