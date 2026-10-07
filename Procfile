release: python backend/manage.py migrate --noinput && python backend/manage.py invalidate_api_cache
worker: celery --workdir backend --app fpbase worker --concurrency 1 --max-memory-per-child 400000 --pool=prefork --without-gossip --without-mingle --without-heartbeat --loglevel=warning
web: gunicorn --chdir backend -c backend/config/gunicorn.py --workers 2 --worker-class sync --max-requests 500 --max-requests-jitter 50 --timeout 30 --access-logfile /dev/null --log-level warning config.wsgi:application
