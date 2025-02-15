release: python backend/manage.py migrate --noinput
worker: celery --workdir backend --app fpbase worker --concurrency 4 --without-gossip --without-mingle --without-heartbeat --loglevel=info
web: gunicorn --chdir backend config.wsgi:application
