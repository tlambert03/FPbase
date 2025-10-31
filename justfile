# list all available just commands
@_default:
    @just --list

# install dependencies and set up the development environment
setup:
    pnpm install
    uv sync

# start both frontend and backend servers
serve:
    pnpm dev

shell:
    uv run backend/manage.py shell_plus

# clone the production database to local
pgpull:
    dropdb fpbase
    heroku  pg:pull DATABASE_URL fpbase --exclude-table-data='public.django_session;public.proteins_ocfluoreff' -a fpbase || true
    uv run backend/manage.py migrate

# start the frontend server
frontend:
    pnpm --stream -r start

# start the backend server
backend:
    uv run backend/manage.py runserver

# run production build locally (build frontend, collect static, run gunicorn with whitenoise)
# uses config.settings.production_local which has all dummy credentials for local testing
prod-local:
    #!/usr/bin/env bash
    set -euo pipefail
    echo "🏗️  Building frontend assets..."
    pnpm build
    echo "📦 Collecting static files..."
    DJANGO_SETTINGS_MODULE=config.settings.production_local \
        uv run backend/manage.py collectstatic --noinput --clear
    echo "🚀 Starting gunicorn on http://localhost:8000"
    echo "   Press Ctrl+C to stop"
    echo "⚠️  Using production-like config with local resources (see config/settings/production_local.py)"
    DJANGO_SETTINGS_MODULE=config.settings.production_local \
        uv run gunicorn --chdir backend --workers 2 --worker-class sync --timeout 30 --bind 127.0.0.1:8000 config.wsgi:application

test-js:
    pnpm --filter @fpbase/spectra test:ci

test-e2e:
    uv run pytest backend/tests_e2e/ -v -n=6 --browser chromium
    uv run pytest backend/tests_e2e/ -v -n=6 --browser webkit

test-py:
    uv run pytest -v -n=6
    uv run pytest backend/tests_e2e/ -v -n=6

snapshots-update:
    uv run pytest backend/tests_e2e/ --visual-snapshots -n 4 --update-snapshots

snapshots-test:
    uv run pytest backend/tests_e2e/ --visual-snapshots -n 4

test: test-py test-js test-e2e

# clean up all virtual environments, caches, and build artifacts
clean:
    find . -name __pycache__ -type d -exec rm -r {} +
    find . -name '*.pyc' -type f -delete
    rm -rf backend/staticfiles
    rm -rf frontend/dist
    rm -rf frontend/node_modules
    rm -rf node_modules
    rm -rf .venv
    rm -rf .pytest_cache
    rm -rf .ruff_cache
    rm -rf .mypy_cache
    rm -f coverage.xml
    rm -rf __snapshots__
    rm -rf snapshot_failures
