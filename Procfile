release: python manage.py migrate --noinput
web: daphne config.asgi:application --port $PORT --bind 0.0.0.0
