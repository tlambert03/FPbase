# Overview

This repo is a Django web app for <https://www.fpbase.org> with:

- PostgreSQL database
- both REST + GraphQL APIs (preferring GraphQL for new code)
- Celery background tasks (using Redis as a broker in production)
- React frontend (but lots of legacy code and a mishmash of patterns being used)
- deployed on Heroku

## Main tooling

- `uv` for all python dependencies (defined in `pyproject.toml`)
- `pnpm` for all JS dependencies (defined in `package.json`, with additional `package.json` files in frontend/ and packages/*)
- `vite` as the frontend build tool (configured in `frontend/vite.config.ts`)
- `biome` for JS/TS linting and formatting (configured in `biome.json`)
- `prek` (modern replacement of pre-commit) for git hooks (configured in `.pre-commit-config.yaml`).
- `pytest` for Python testing, which covers both backend and frontend (end-to-end tests using Playwright)
- `just` for task running (configured in `justfile`)
- `django_vite` is used to add our Vite-built frontend assets into Django templates
- `gh` CLI for querying the GitHub API from the command line (prefer this over MCP)

## Bash commands

- `pnpm install && uv sync`: Install all dependencies - call whenever you checkout a new branch or pull changes
- `uv run pytest`: Run Python backend unit tests
- `uv run pytest backend/tests_e2e/ -n 6`: Run end-to-end tests with playwright.
  *this will autobuild the frontend first if needed!  no need to call pnpm build first*
- `pnpm build`: Build static frontend assets
- `pnpm dev`: Run *both* the React frontend dev server with HMR and Django backend dev server.
  *use this for rapid development! and visit localhost:8000 with playwright. vite takes care of rebuilding assets as needed*

- `prek -a`: Run all pre-commit hooks on all files (takes care of linting/formatting for both Python and JS/TS)
- `uv run basedpyright`: Run Python typecheck
- `pnpm typecheck`: Run TypeScript typecheck

- `uv run backend/manage.py shell_plus` : Open a Django shell with all models pre-imported

## Code style

- Use type hints in Python code wherever possible (though much of the codebase is not yet typed)
- prefer functional pytest style over class-based tests
- Use ES modules (import/export) syntax, not CommonJS (require)
- Destructure imports when possible (eg. import { foo } from 'bar')
- Use biome for JS/TS linting and formatting
- Use ruff for Python linting and formatting
- DO NOT use `time.sleep()` or `page.wait_for_timeout()` in tests; use proper waits on a specific condition instead.
- avoid nested imports unless specifically used to avoid circular imports or delay heavy imports.

## Workflow

- NEVER commit directly to main branch
- Be sure to typecheck when you’re done making a series of code changes
- Be sure to run all tests (unit + e2e) before pushing code
- Tests use `--reuse-db` by default (from pyproject.toml). If you see DB errors in tests, try rerunning, or `uv run pytest --create-db`.
- When fixing Sentry issues, do NOT just silence errors - you MUST fix the root cause!
- the current year is 2025 (not 2024), for web-searches

## Icon System

FPbase uses an inline SVG icon system with library abstraction:
- Icons are rendered via `{% icon "name" %}` template tag
- SVG data is extracted from npm packages (FontAwesome, Lucide) into JSON files
- To switch libraries: set `FPBASE_ICON_LIBRARY=lucide` environment variable
- To add new icons: update `scripts/extract_icons.py` mappings and run `python scripts/extract_icons.py all`
- See `docs/icons.md` for full documentation

## Agent tips

- If you don't know something, prioritize using WebSearch (or the context7 mcp) to find the answer.
  Don't make up answers or waste too much time guessing and checking.
- when writing code, don't be overly defensive unless the application warrants it. don't add
  extra error handling "just in case".
