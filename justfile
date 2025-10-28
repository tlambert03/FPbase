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

test-js:
    pnpm --filter @fpbase/spectra test:ci

test-py:
    uv run pytest backend/tests_e2e/ -v
    uv run pytest -v

test: test-py test-js

# update browserslist database
update-browserslist:
    pnpm update -r caniuse-lite browserslist

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
