release: python manage.py migrate --noinput
web: newrelic-admin run-program daphne config.asgi:application --port $port --bind 0.0.0.0
