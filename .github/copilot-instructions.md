# FPbase Copilot Instructions

## Project Overview
FPbase is a Django/React monorepo for the Fluorescent Protein Database (fpbase.org). The backend uses Django 4.x with GraphQL and REST APIs, while the frontend is a hybrid of server-rendered Django templates with embedded React apps built via Webpack and Vite.

## Architecture

### Monorepo Structure
- **`backend/`**: Django application (Python 3.11+)
  - `fpbase/`: Core Django app with users, views, and utilities
  - `proteins/`: Main domain logic for proteins, spectra, microscopes, and optical configurations
  - `references/`: Publication and citation management
  - `fpseq/`: Bioinformatics sequence alignment and analysis using Biopython
  - `favit/`: User favorites system
- **`frontend/`**: Webpack-bundled assets integrated via `django-webpack-loader`
- **`packages/`**: Standalone Vite apps (`blast/`, `spectra/`) embedded in Django templates

### Key Technologies
- **Backend**: Django 4.x, Django REST Framework, Graphene (GraphQL), Celery (Redis), PostgreSQL
- **Frontend**: React 16, Material-UI, Webpack 5, Vite 4
- **Search**: Algolia for protein/organism search
- **Bioinformatics**: Biopython (sequence alignment), BLAST (local binaries in `backend/bin/`)
- **Deployment**: Heroku, AWS S3 (media), Sentry (error tracking)

## Development Workflow

### Setup (from README.md)
```bash
# Python 3.13 environment
uv sync
# PostgreSQL setup
createdb fpbase
python backend/manage.py migrate
# Frontend
pnpm install
# Start dev servers
npm run start  # or: pnpm dev (runs both frontend + backend)
python backend/manage.py runserver
```

### Running Tests
```bash
# Backend tests (pytest)
uv run pytest --color=yes -v --cov
# Tests require PostgreSQL and Redis (see .github/workflows/ci.yml)
```

### Settings Architecture
Django settings split across `backend/config/settings/`:
- `base.py`: Shared configuration
- `local.py`: Development (DEBUG=True, dummy cache, console email)
- `production.py`: Production (Heroku, AWS S3, real cache)
- `test.py`: Testing (uses `MockWebpackLoader` to skip frontend builds)

Environment variables loaded via `django-environ` from `.env` file (set `DJANGO_READ_DOT_ENV_FILE=True`).

## Domain Model Patterns

### Proteins App (`backend/proteins/`)
Core models in `models/`:
- `Protein`: The main FP entity with sequence data, optical properties, lineage tracking
- `State`/`Fluorophore`/`Dye`: Different fluorescent entities with associated `Spectrum` data
- `Microscope`/`OpticalConfig`: Microscope configurations with filters and calculated efficiencies
- `Spectrum`: Spectral data (excitation/emission) with CSV data storage
- Use `mixins.Authorable` for models with creator/owner tracking

### Key Conventions
1. **Sequence Analysis**: Use `fpseq.FPSeq` class (wrapper around Biopython) for protein sequences
   - Alignment via `fpseq.align.align_seqs()` using Biopython
   - Mutations tracked via `fpseq.mutations.get_mutations()`
2. **Factories**: Use `factory_boy` for test fixtures (see `proteins/factories.py`, `references/factories.py`)
   - Example: `ProteinFactory`, `SpectrumFactory`, `OrganismFactory`
3. **Slugs**: All public-facing models use auto-generated slugs for URLs
4. **Versioning**: `django-reversion` tracks changes to proteins/microscopes

### Search & Indexing
- Algolia integration via `algoliasearch_django` (see `proteins/index.py`)
- Protein and Organism models auto-indexed on save
- Use `ALGOLIA_SUFFIX` env var to separate dev/prod indexes

## API Patterns

### REST API (`/api/`)
- DRF viewsets in `proteins/api/` with `drf-spectacular` for OpenAPI docs
- URL config in `backend/config/api_router.py`
- Use `djangorestframework-csv` for CSV export endpoints

### GraphQL (`/graphql/`)
- Schema defined in `fpbase/schema.py` (imports from `proteins.schema` and `references.schema`)
- Uses `graphene-django` with `graphene-django-optimizer` for query optimization
- Auth via `django-graphql-jwt`

## Frontend Integration

### Webpack + Django Templates
- Webpack builds to `frontend/dist/` with stats tracked in `webpack-stats.json`
- Django loads bundles via `{% load webpack_loader %}` template tags
- Hot reload available with `HOT_RELOAD=1` env var
- Entry points in `frontend/src/`: `index.js`, `spectra-viewer.js`, `blast-app.js`, etc.

### Vite Apps (packages/)
- Standalone React apps (`@fpbase/blast`, `@fpbase/spectra`) embedded as iframes or via CDN
- Built separately with `pnpm -r build`

## Background Tasks
- Celery tasks in `proteins/tasks.py` for long-running calculations
  - `calc_fret()`: FRET efficiency calculations
  - `calculate_scope_report()`: Microscope optical efficiency reports
- Redis broker required for local dev

## Testing Patterns
- Tests in `*/tests/` directories (pytest)
- Use `@pytest.mark.django_db` for database access
- Frontend-dependent tests use `@pytest.mark.usefixtures("uses_frontend", "use_real_webpack_loader")`
- Factory fixtures preferred over manual object creation

## Custom Middleware
- `fpbase.middleware.CanonicalDomainMiddleware`: Redirects to canonical domain in production
- `fpbase.middleware.BlackListMiddleware`: IP blocking (uses `BLOCKED_IPS` setting)

## Common Commands
```bash
# Migrations
python backend/manage.py makemigrations
python backend/manage.py migrate

# Create superuser
python backend/manage.py createsuperuser

# Frontend build
pnpm build  # production build
pnpm start  # dev server with watch

# Algolia reindex
python backend/manage.py algolia_reindex

# Django shell
python backend/manage.py shell_plus  # (django-extensions)
```

## Code Quality
- Uses `ruff` for Python linting/formatting (config in `pyproject.toml`)
- Pre-commit hooks configured (see `.pre-commit-config.yaml` if exists)
- Coverage tracked via `pytest-cov` and codecov.io

## Important Files
- `backend/config/urls.py`: Main URL routing
- `backend/fpbase/sitemaps.py`: SEO sitemaps for proteins, organisms, microscopes
- `backend/proteins/util/helpers.py`: Utility functions (slug caching, spectra plotting, FRET calculations)
- `pyproject.toml`: Python dependencies and pytest config
- `pnpm-workspace.yaml`: Monorepo workspace config
