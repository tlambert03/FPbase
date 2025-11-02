# FPbase - Fluorescent Protein Database

Django web app for <https://www.fpbase.org> with React frontend. PostgreSQL database, REST + GraphQL APIs, Celery background tasks.

## Important instructions

- don't auto-deploy to heroku from this local git repo.  All deployments must go through github PRs.
- unless explicitly told otherwise, don't run `git commit`, let me review changes.
- avoid nested imports unless specifically used to avoid circular imports or delay heavy imports.
- avoid arbitrary time.sleep() calls in tests; use proper waits on a specific condition instead.

## Tech Stack

**Backend**: Django, Python, DRF, GraphQL (graphene-django), PostgreSQL, Celery + Redis
**Frontend**: React, Webpack, pnpm monorepo (packages: spectra, blast use Vite)
**Science**: BioPython, NumPy, Pandas, SciPy, Matplotlib

## Key Overrides

- **Line length: 119 chars** (not 88) - configured in pyproject.toml
- Frontend uses pnpm, not npm
- Python uses `uv` for virtualenv management and running commands

## Common Commands

```bash
# Setup
uv sync                                # Install/update Python deps
pnpm install                           # Install Node deps

# Development
pnpm dev                               # Start webpack + Django dev server
uv run backend/manage.py shell_plus    # Interactive shell with auto-imports

# Testing
uv run pytest path/to/test.py          # Run specific test
uv run pytest --cov                    # Run with coverage
pnpm --filter @fpbase/spectra test:ci  # Frontend tests

# Code Quality
ruff format backend                    # Format Python
ruff check --fix backend               # Auto-fix linting
uv run mypy backend                    # Type check
prek run --all-files                   # Run all hooks

# Build
pnpm build                             # Build all frontend packages
```

## Testing Notes

- Tests colocated with apps: `proteins/tests/`, `fpbase/tests/`
- pytest with Django plugin, --reuse-db enabled
- factory-boy for test data
- Selenium + Playwright for browser tests
- Settings: `config.settings.test`
- Coverage target: 100% (soft), PRs: 5% for new code

## Important Tools

- **django-extensions**: `shell_plus` auto-imports models
- **django-reversion**: Model history tracking enabled
- **Algolia**: Full-text search integration
- **Sentry**: Error tracking (production)
- **Celery**: Background tasks (e.g., BLAST searches)

## Development Notes

- Django settings via django-environ (loads from .env)
- Pre-commit hooks auto-run (ruff format, ruff check --fix, django-upgrade)
- Database: PostgreSQL required (not SQLite)
- Dual APIs: REST (drf-spectacular docs) + GraphQL (JWT auth)
- Static files: WhiteNoise with Brotli compression
- Deployment: Heroku (Procfile: web=gunicorn, worker=celery)

## Gotchas

- Migrations in proteins/ are extensive - review carefully before changing models
