release: pdm install && pdm run migrate --noinput
worker: celery --app fpbase worker --concurrency 4 --without-gossip --without-mingle --without-heartbeat
web: pdm run gunicorn config.wsgi:application
