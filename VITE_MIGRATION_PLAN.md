# Webpack to Vite 7 Migration Plan for FPbase

## GENERAL INSTRUCTIONS FOR CLAUDE CODE

### Migration Overview

Migrate FPbase's frontend build system from Webpack 5 + Babel to Vite 7, leveraging django-vite for Django integration.

### Decision Summary

- ✅ **Migrate to Vite 7:** YES (Vite is the modern standard, significant performance gains)
- ✅ **Keep 7 Entry Points:** YES (Required for D3 v3/v7 isolation)
- ✅ **Template Migration:** Big Bang (All at once, cleaner approach)

### Tips and instructions

- work step-by-step, and when you are happy with a specific step, commit it to the `use-vite` branch.
- do NOT push to remote, just commit locally.
- run tests frequently with `just test` (which will test both backend and frontend).
- if you need to reinstall dependencies, run `just setup`. (which doed `pnpm install` and `uv sync`)
- if you need to start a dev server, run `pnpm dev` from the root of the monorepo (it will start both vite and django dev servers).
- if you want to interact the dev site in a browser, use playwright mcp and open `http://localhost:8000`.
- prefer RAG and websearch over guess and check
- refer frequently to the official documentation with questions:
  - Vite: <https://vite.dev/guide/> and <https://vite.dev/config/>
  - django-vite: <https://github.com/MrBin99/django-vite>

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current Architecture Analysis](#current-architecture-analysis)
3. [Target Architecture](#target-architecture)
4. [Migration Strategy](#migration-strategy)
5. [Detailed Implementation Plan](#detailed-implementation-plan)
6. [Risk Assessment & Mitigation](#risk-assessment--mitigation)
7. [Testing Strategy](#testing-strategy)
8. [Reference Information](#reference-information)

---

## Current Architecture Analysis

### 1. Webpack Configuration

**File:** `frontend/webpack.config.js`

#### Entry Points (7 total)

```javascript
entry: {
  main: './src/index.js',              // Core site functionality
  embedscope: './src/embedscope.js',   // Microscope embed (uses D3 v3 CDN)
  litemol: './src/my-litemol.js',      // Protein structure viewer
  spectraViewer: './src/spectra-viewer.js',       // Full spectra viewer
  simpleSpectraViewer: './src/simple-spectra-viewer.js',  // Simple spectra
  microscopeForm: './src/microscope-form.js',     // Microscope config form
  blast: './src/blast-app.js',         // BLAST search
  proteinTable: './src/protein-table.js',  // Protein table
}
```

#### Code Splitting Strategy

**6-tier manual chunking:**

```javascript
splitChunks: {
  cacheGroups: {
    sentry: { priority: 30 },        // Shared error tracking
    react: { priority: 25 },         // React + ReactDOM + scheduler
    jquery: { priority: 20 },        // jQuery (used across bundles)
    d3: { priority: 15 },            // D3 v7 (NOT in embedscope)
    vendors: { priority: 10 },       // Other node_modules
    common: { priority: 5 },         // Shared app code
  }
}
```

**CRITICAL:** `embedscope` must NOT include D3 v7 (uses CDN D3 v3 instead).

#### Asset Processing

- **Sass/SCSS** → CSS with autoprefixer + cssnano
- **JavaScript** → Babel (preset-env, preset-react)
- **Assets** → Copy plugin for `microscope.js`
- **Source Maps** → Sentry webpack plugin uploads

#### Dev Server

- Port: 8080
- Hot Reload: `HOT_RELOAD=1` env var
- CORS enabled for Django backend
- Public path: `http://localhost:8080/static/`

### 2. Django Integration

#### Settings Configuration

**File:** `backend/config/settings/base.py` (lines 226-237)

```python
INSTALLED_APPS.append("webpack_loader")

WEBPACK_LOADER = {
    "DEFAULT": {
        "CACHE": not DEBUG,
        "BUNDLE_DIR_NAME": "/",
        "STATS_FILE": str(ROOT_DIR.parent / "frontend" / "dist" / "webpack-stats.json"),
        "POLL_INTERVAL": 0.1,
        "TIMEOUT": None,
        "IGNORE": [r".*\.hot-update.js", r".+\.map"],
    }
}
```

**File:** `backend/config/settings/test.py` (lines 100-105)

```python
class MockWebpackLoader(FakeWebpackLoader):
    def get_assets(self):
        return {}

WEBPACK_LOADER["DEFAULT"]["LOADER_CLASS"] = "config.settings.test.MockWebpackLoader"
```

#### Template Usage

**10 templates use `render_bundle`:**

| Template | Bundles Used | Notes |
|----------|-------------|-------|
| `fpbase/templates/base.html` | main (css + js) | Base template for all pages |
| `proteins/templates/spectra.html` | spectraViewer | Full spectra viewer page |
| `proteins/templates/spectra_graph.html` | simpleSpectraViewer | Embedded spectra widget |
| `proteins/templates/table.html` | proteinTable | Protein comparison table |
| `proteins/templates/compare.html` | simpleSpectraViewer | Protein comparison |
| `proteins/templates/proteins/microscope_form.html` | microscopeForm | Microscope config |
| `proteins/templates/proteins/microscope_embed.html` | main (css) + embedscope (js) | **CRITICAL: D3 v3** |
| `proteins/templates/proteins/protein_detail.html` | simpleSpectraViewer + litemol | Protein page |
| `proteins/templates/proteins/blast.html` | blast | BLAST search |
| `fpbase/templates/500.html` | main (css) | Error page |

**Pattern Examples:**

```django
{% load render_bundle from webpack_loader %}
{% render_bundle 'main' 'css' %}
{% render_bundle 'main' 'js' %}
{% render_bundle 'embedscope' 'js' %}
{% render_bundle 'litemol' attrs='defer' %}
```

### 3. E2E Test Infrastructure

**File:** `backend/tests_e2e/conftest.py` (lines 121-142)

```python
@pytest.fixture(scope="module", autouse=True)
def _setup_frontend_assets() -> None:
    """Build webpack assets once per test module if needed."""
    stats_file = Path(django.conf.settings.WEBPACK_LOADER["DEFAULT"]["STATS_FILE"])

    if _frontend_assets_need_rebuild(stats_file):
        print("Building frontend assets for e2e tests...")
        subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)

    # UNDO MockWebpackLoader for e2e tests
    utils.get_loader = _get_real_get_loader
```

**Logic:**

- Checks if `webpack-stats.json` exists and is valid
- Checks if source files are newer than stats file
- Rebuilds frontend if needed
- Reverts MockWebpackLoader to real loader

### 4. Package Structure

**Root monorepo:** pnpm workspaces

```
fpbase/
├── frontend/              # Main webpack bundle (this migration)
│   ├── dist/              # Build output
│   ├── src/
│   │   ├── index.js       # main entry
│   │   ├── embedscope.js  # embedscope entry (D3 v3 CDN)
│   │   ├── blast-app.js   # blast entry
│   │   └── ...
│   ├── package.json
│   └── webpack.config.js
├── packages/
│   ├── spectra/           # ✅ Already uses Vite 4
│   │   ├── vite.config.js
│   │   └── package.json
│   ├── blast/             # ✅ Already uses Vite 4
│   │   ├── vite.config.js
│   │   └── package.json
│   └── protein-table/     # ✅ Already uses Vite 4
└── package.json
```

**Note:** `spectra`, `blast`, and `protein-table` are consumed by webpack as workspace packages. In Vite, they'll be directly importable.

### 5. Key Dependencies

**From `frontend/package.json`:**

```json
{
  "dependencies": {
    "@fpbase/blast": "workspace:*",
    "@fpbase/spectra": "workspace:*",
    "@sentry/browser": "^10.22.0",
    "jquery": "^3.7.0",
    "d3": "^7.9.0",
    "react": "^19.2.0",
    "react-dom": "^19.2.0",
    // ... bootstrap, select2, etc.
  },
  "devDependencies": {
    "webpack": "^5.102.1",
    "webpack-bundle-tracker": "3.2.1",
    "babel-loader": "^10.0.0",
    // ... all webpack/babel plugins
  }
}
```

---

## Target Architecture

### 1. Vite Configuration

**New file:** `frontend/vite.config.js`

```javascript
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { viteStaticCopy } from 'vite-plugin-static-copy';
import { sentryVitePlugin } from '@sentry/vite-plugin';

export default defineConfig(({ mode }) => {
  const isDev = mode === 'development';

  return {
    // Match Django STATIC_URL
    base: '/static/',

    // Build configuration
    build: {
      // Output to frontend/dist/
      outDir: 'dist',
      emptyOutDir: true,

      // Generate manifest.json for django-vite
      manifest: 'manifest.json',

      // Source maps for Sentry
      sourcemap: true,

      // Multi-page app configuration
      rollupOptions: {
        input: {
          main: './src/index.js',
          embedscope: './src/embedscope.js',
          litemol: './src/my-litemol.js',
          spectraViewer: './src/spectra-viewer.js',
          simpleSpectraViewer: './src/simple-spectra-viewer.js',
          microscopeForm: './src/microscope-form.js',
          blast: './src/blast-app.js',
          proteinTable: './src/protein-table.js',
        },

        // Manual code splitting (hybrid approach)
        output: {
          manualChunks(id) {
            // Sentry (shared across all)
            if (id.includes('node_modules/@sentry')) {
              return 'vendor-sentry';
            }

            // React (shared, but NOT in embedscope)
            if (id.includes('node_modules/react')) {
              return 'vendor-react';
            }

            // jQuery (shared across multiple)
            if (id.includes('node_modules/jquery')) {
              return 'vendor-jquery';
            }

            // D3 v7 (EXCLUDE from embedscope - it uses CDN D3 v3)
            if (id.includes('node_modules/d3')) {
              // Check if this is being imported by embedscope
              // If so, throw error to prevent bundling
              return 'vendor-d3';
            }
          },
        },
      },
    },

    // Development server
    server: {
      port: 5173, // Vite default (change if needed)
      strictPort: true,
      origin: 'http://localhost:5173',
      cors: true,
      hmr: {
        protocol: 'ws',
        host: 'localhost',
      },
    },

    // Resolve aliases (match webpack)
    resolve: {
      alias: {
        '@fpbase/spectra': '/packages/spectra/src/index.jsx',
        '@fpbase/blast': '/packages/blast/src/index.js',
        '@fpbase/protein-table': '/packages/protein-table/src/index.jsx',
      },
    },

    // Plugins
    plugins: [
      // React with Fast Refresh
      react(),

      // Copy static files (microscope.js)
      viteStaticCopy({
        targets: [
          {
            src: '../backend/fpbase/static/js/microscope.js',
            dest: 'js',
          },
        ],
      }),

      // Sentry source map upload (production only)
      !isDev && sentryVitePlugin({
        org: 'talley-lambert',
        project: 'fpbase',
        authToken: process.env.SENTRY_AUTH_TOKEN,
        release: process.env.HEROKU_SLUG_COMMIT,
      }),
    ].filter(Boolean),

    // CSS configuration
    css: {
      postcss: {
        plugins: [
          require('autoprefixer'),
          require('cssnano'),
        ],
      },
    },

    // Define environment variables
    define: {
      'process.env.NODE_ENV': JSON.stringify(mode),
      'process.env.SENTRY_DSN': JSON.stringify(process.env.SENTRY_DSN || ''),
      'process.env.HEROKU_SLUG_COMMIT': JSON.stringify(process.env.HEROKU_SLUG_COMMIT || ''),
    },
  };
});
```

### 2. Django Settings Updates

**File:** `backend/config/settings/base.py`

```python
# Remove webpack_loader, add django_vite
INSTALLED_APPS = [
    'django_vite',  # Add BEFORE other apps
    # ... rest of INSTALLED_APPS
]

# Remove WEBPACK_LOADER, add DJANGO_VITE
DJANGO_VITE = {
    "default": {
        "dev_mode": DEBUG,  # Match DEBUG setting
        "dev_server_protocol": "http",
        "dev_server_host": "localhost",
        "dev_server_port": 5173,
        "static_url_prefix": "",
        "manifest_path": str(ROOT_DIR.parent / "frontend" / "dist" / "manifest.json"),
    }
}
```

**File:** `backend/config/settings/test.py`

```python
# Remove MockWebpackLoader

# django-vite in test mode uses manifest
DJANGO_VITE = {
    "default": {
        "dev_mode": False,  # Always use manifest in tests
        "manifest_path": str(ROOT_DIR.parent / "frontend" / "dist" / "manifest.json"),
    }
}
```

### 3. Template Migration

**Pattern changes:**

| Old (webpack_loader) | New (django_vite) |
|---------------------|-------------------|
| `{% load render_bundle from webpack_loader %}` | `{% load django_vite %}` |
| `{% render_bundle 'main' 'css' %}` | `{% vite_asset 'src/index.js' %}` (auto-includes CSS) |
| `{% render_bundle 'main' 'js' %}` | `{% vite_asset 'src/index.js' %}` |
| `{% render_bundle 'embedscope' 'js' %}` | `{% vite_asset 'src/embedscope.js' %}` |
| `{% render_bundle 'litemol' attrs='defer' %}` | `{% vite_asset 'src/my-litemol.js' defer %}` |

**New addition (for HMR in dev):**

Add to `<head>` in base.html:

```django
{% load django_vite %}
{% vite_hmr_client %}  <!-- Only in dev mode -->
```

**Complete migration list:**

1. `backend/fpbase/templates/base.html`
   - Add `{% vite_hmr_client %}` in `<head>`
   - Replace `{% render_bundle 'main' 'css' %}` → `{% vite_asset 'src/index.js' %}`
   - Replace `{% render_bundle 'main' 'js' %}` → (already included above)

2. `backend/proteins/templates/spectra.html`
   - Replace `{% render_bundle 'spectraViewer' %}` → `{% vite_asset 'src/spectra-viewer.js' %}`

3. `backend/proteins/templates/spectra_graph.html`
   - Replace `{% render_bundle 'simpleSpectraViewer' %}` → `{% vite_asset 'src/simple-spectra-viewer.js' %}`

4. `backend/proteins/templates/table.html`
   - Replace `{% render_bundle 'proteinTable' %}` → `{% vite_asset 'src/protein-table.js' %}`

5. `backend/proteins/templates/compare.html`
   - Replace `{% render_bundle 'simpleSpectraViewer' %}` → `{% vite_asset 'src/simple-spectra-viewer.js' %}`

6. `backend/proteins/templates/proteins/microscope_form.html`
   - Replace `{% render_bundle 'microscopeForm' %}` → `{% vite_asset 'src/microscope-form.js' %}`

7. `backend/proteins/templates/proteins/microscope_embed.html` **[CRITICAL: D3 v3]**
   - Replace `{% render_bundle 'main' 'css' %}` → `{% vite_asset 'src/index.js' %}`
   - Replace `{% render_bundle 'embedscope' 'js' %}` → `{% vite_asset 'src/embedscope.js' %}`
   - **VERIFY:** embedscope.js must NOT import D3 v7

8. `backend/proteins/templates/proteins/protein_detail.html`
   - Replace `{% render_bundle 'simpleSpectraViewer' %}` → `{% vite_asset 'src/simple-spectra-viewer.js' %}`
   - Replace `{% render_bundle 'litemol' attrs='defer' %}` → `{% vite_asset 'src/my-litemol.js' defer %}`

9. `backend/proteins/templates/proteins/blast.html`
   - Replace `{% render_bundle 'blast' %}` → `{% vite_asset 'src/blast-app.js' %}`

10. `backend/fpbase/templates/500.html`
    - Replace `{% render_bundle 'main' 'css' %}` → `{% vite_asset 'src/index.js' %}`

### 4. E2E Test Updates

**File:** `backend/tests_e2e/conftest.py`

```python
@pytest.fixture(scope="module", autouse=True)
def _setup_frontend_assets() -> None:
    """Build Vite assets once per test module if needed."""
    # Change from webpack-stats.json to manifest.json
    manifest_file = Path(django.conf.settings.DJANGO_VITE["default"]["manifest_path"])

    if _frontend_assets_need_rebuild(manifest_file):
        print("Building frontend assets for e2e tests...")
        # Change from webpack build to vite build
        subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)

    # No need to revert MockWebpackLoader (django-vite doesn't need mocking)
```

Update `_frontend_assets_need_rebuild()`:

```python
def _frontend_assets_need_rebuild(manifest_file) -> bool:
    """Check if frontend assets need to be rebuilt."""
    if not manifest_file.is_file():
        return True

    # Check if manifest is valid JSON
    try:
        manifest = json.loads(manifest_file.read_bytes())
        if not manifest:
            return True
    except (json.JSONDecodeError, ValueError):
        return True

    # Check if source files are newer
    manifest_mtime = manifest_file.stat().st_mtime
    frontend_src = Path(__file__).parent.parent.parent / "frontend" / "src"
    if any(
        f.stat().st_mtime > manifest_mtime for f in frontend_src.rglob("*")
        if f.is_file() and not f.name.startswith(".")
    ):
        return True

    return False
```

---

## Migration Strategy

### Phase Overview

```
Phase 1: Preparation
├─ Install dependencies
├─ Create vite.config.js
├─ Test standalone build
└─ Verify manifest.json

Phase 2: Vite Configuration
├─ Configure all 7 entry points
├─ Verify D3 v3/v7 isolation
├─ Test dev server + HMR
└─ Profile bundle sizes

Phase 3: Django Integration
├─ Update settings (base.py, test.py)
├─ Migrate 10 templates
├─ Test each page loads
└─ Verify no console errors

Phase 4: Testing & Validation
├─ Update E2E test setup
├─ Run full test suite
├─ Manual QA on staging
└─ Fix any issues

Phase 5: Production Deployment
├─ Remove webpack dependencies
├─ Update documentation
├─ Deploy to production
└─ Monitor for issues

```

## Detailed Implementation Plan

### Phase 1: Preparation

#### Step 1.1: Install Vite Dependencies

```bash
cd frontend

# Remove webpack dependencies (will be done in Phase 5)
# For now, keep both to allow rollback

# Install Vite and plugins
pnpm add -D vite@^7.1.12 \
  @vitejs/plugin-react@^5.1.0 \
  vite-plugin-static-copy@^3.1.4

# Install Sentry Vite plugin
pnpm add -D @sentry/vite-plugin@^4.6.0

# Verify Node.js version (Vite 7 requires Node 20.19+ or 22.12+)
node --version  # Should be 20.19+ or 22.12+
```

#### Step 1.2: Create vite.config.js

Create `frontend/vite.config.js` with the configuration from [Target Architecture](#1-vite-configuration).

#### Step 1.3: Update package.json Scripts

```json
{
  "scripts": {
    "dev": "vite",
    "dev:old": "HOT_RELOAD=1 webpack-dev-server --mode development",
    "build": "vite build",
    "build:old": "NODE_ENV=production webpack --mode production",
    "preview": "vite preview"
  }
}
```

#### Step 1.4: Add Modulepreload Polyfill

**Edit `frontend/src/index.js`** (and other entry points):

```javascript
// Add as first import
import 'vite/modulepreload-polyfill';

// Rest of imports...
```

Do this for all 8 entry points:

- `src/index.js`
- `src/embedscope.js`
- `src/blast-app.js`
- `src/spectra-viewer.js`
- `src/simple-spectra-viewer.js`
- `src/microscope-form.js`
- `src/protein-table.js`
- `src/my-litemol.js`

#### Step 1.5: Test Standalone Vite Build

```bash
cd frontend

# Clean old build
rm -rf dist

# Run Vite build
pnpm run build

# Verify output
ls -lah dist/

# Expected files:
# - manifest.json
# - assets/main-[hash].js
# - assets/main-[hash].css
# - assets/embedscope-[hash].js
# - ... (all entry points)
# - assets/vendor-sentry-[hash].js
# - assets/vendor-react-[hash].js
# - assets/vendor-jquery-[hash].js
# - assets/vendor-d3-[hash].js
# - js/microscope.js (static copy)

# Verify manifest.json structure
cat dist/manifest.json | jq .
```

**Success Criteria:**

- ✅ Build completes without errors
- ✅ `manifest.json` exists and is valid JSON
- ✅ All 8 entry points have generated JS files
- ✅ CSS files are extracted
- ✅ Vendor chunks are created
- ✅ `microscope.js` is copied to `js/`

---

### Phase 2: Vite Configuration

#### Step 2.1: Verify D3 v3/v7 Isolation

**CRITICAL:** `embedscope.js` uses CDN D3 v3, must NOT bundle D3 v7.

**Check embedscope imports:**

```bash
cd frontend
grep -r "from 'd3'" src/embedscope.js
grep -r "import.*d3" src/embedscope.js
```

**Expected:** NO matches. Embedscope should not import D3.

**If D3 is imported:**

1. Refactor embedscope to remove D3 v7 dependency
2. Ensure it only uses `window.d3` (CDN D3 v3)

#### Step 2.2: Profile Bundle Sizes

```bash
cd frontend

# Build with analysis
pnpm run build

# Analyze output
du -sh dist/assets/*.js | sort -h

# Compare to webpack baseline
# Expected: ±10% variance is acceptable
# If >20% larger, investigate and optimize
```

**Baseline (webpack):** (record actual sizes from current build)

```
main: ~XXX KB
embedscope: ~XXX KB
blast: ~XXX KB
...
```

#### Step 2.3: Test Dev Server + HMR

```bash
cd frontend

# Start Vite dev server
pnpm run dev

# Server should start on http://localhost:5173
# Open browser to http://localhost:5173/src/index.js
# (Will NOT work yet - Django integration needed)

# Test HMR by editing a file
# Edit src/index.js, save
# Check terminal for HMR update message
```

**Success Criteria:**

- ✅ Dev server starts successfully
- ✅ HMR updates reflect in browser (<1s)
- ✅ No errors in browser console
- ✅ No errors in terminal

#### Step 2.4: Test Sentry Integration

**Verify source maps:**

```bash
cd frontend

# Build
pnpm run build

# Check for .map files
ls -lah dist/assets/*.map

# Verify Sentry plugin ran (if SENTRY_AUTH_TOKEN set)
# Check build output for "Sentry: Uploading source maps"
```

---

### Phase 3: Django Integration

#### Step 3.1: Install django-vite

```bash
cd /Users/talley/dev/self/FPbase
uv add django-vite
```

#### Step 3.2: Update Django Settings

**Edit `backend/config/settings/base.py`:**

```python
# Line 226: Change from webpack_loader to django_vite
INSTALLED_APPS = [
    # ... existing apps
    'django_vite',  # ADD THIS (before apps that use it)
    # 'webpack_loader',  # REMOVE THIS
]

# Lines 228-237: Replace WEBPACK_LOADER with DJANGO_VITE
DJANGO_VITE = {
    "default": {
        "dev_mode": DEBUG,
        "dev_server_protocol": "http",
        "dev_server_host": "localhost",
        "dev_server_port": 5173,
        "static_url_prefix": "",
        "manifest_path": str(ROOT_DIR.parent / "frontend" / "dist" / "manifest.json"),
    }
}

# Remove old WEBPACK_LOADER config
```

**Edit `backend/config/settings/test.py`:**

```python
# Remove lines 10, 100-105 (MockWebpackLoader)

# Add DJANGO_VITE override
DJANGO_VITE = {
    "default": {
        "dev_mode": False,
        "manifest_path": str(ROOT_DIR.parent / "frontend" / "dist" / "manifest.json"),
    }
}
```

#### Step 3.3: Migrate Templates

**Semi-automated approach using sed:**

```bash
cd backend

# Backup templates first
tar -czf templates_backup_$(date +%Y%m%d).tar.gz */templates/

# Step 1: Replace load statement
find . -name "*.html" -exec sed -i '' \
  's/{%\s*load\s*render_bundle\s*from\s*webpack_loader\s*%}/{%\ load\ django_vite\ %}/g' {} \;

# Step 2: Replace render_bundle calls
# Note: These are approximate patterns, verify each manually

# main bundle (most common)
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'main'\s*'css'\s*%}/{%\ vite_asset\ 'src\/index.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'main'\s*'js'\s*%}//g" {} \;  # Remove (auto-included)
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'main'\s*%}/{%\ vite_asset\ 'src\/index.js'\ %}/g" {} \;

# Other bundles
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'embedscope'\s*'js'\s*%}/{%\ vite_asset\ 'src\/embedscope.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'blast'\s*%}/{%\ vite_asset\ 'src\/blast-app.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'spectraViewer'\s*%}/{%\ vite_asset\ 'src\/spectra-viewer.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'simpleSpectraViewer'\s*%}/{%\ vite_asset\ 'src\/simple-spectra-viewer.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'microscopeForm'\s*%}/{%\ vite_asset\ 'src\/microscope-form.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'proteinTable'\s*%}/{%\ vite_asset\ 'src\/protein-table.js'\ %}/g" {} \;
find . -name "*.html" -exec sed -i '' \
  "s/{%\s*render_bundle\s*'litemol'\s*attrs='defer'\s*%}/{%\ vite_asset\ 'src\/my-litemol.js'\ defer\ %}/g" {} \;

# Verify no old tags remain
grep -r "render_bundle" . --include="*.html"
# Expected: No matches (or only in backups/docs)
```

**Manual verification required for each template!**

#### Step 3.4: Add HMR Client to Base Template

**Edit `backend/fpbase/templates/base.html`:**

```django
{% load i18n %}
{% load static %}
{% load django_vite %}  {# Changed from webpack_loader #}
{% load webp_picture from protein_tags %}

<!DOCTYPE html>
<html lang="en">
  <head>
    {% vite_hmr_client %}  {# ADD THIS - enables HMR in dev #}

    {% block ga %}
      ...
    {% endblock ga %}

    ...

    {% vite_asset 'src/index.js' %}  {# Changed from render_bundle 'main' #}

    ...
  </head>
  ...
</html>
```

#### Step 3.5: Test Each Page

For this, you can use the dev server and playwright mcp

**With Django dev server running:**

```bash
# Terminal 1: Both django and vite dev concurrently
pnpm dev

# Open browser to http://localhost:8000
# Test each page:
```

**Checklist:**

- [ ] Homepage loads (main bundle)
- [ ] Protein detail page (simpleSpectraViewer + litemol)
- [ ] Spectra viewer page (spectraViewer)
- [ ] BLAST page (blast)
- [ ] Protein table page (proteinTable)
- [ ] Microscope form page (microscopeForm)
- [ ] Microscope embed page **[CRITICAL: D3 v3]** (embedscope)
- [ ] Compare page (simpleSpectraViewer)
- [ ] Spectra graph widget (simpleSpectraViewer)
- [ ] Error page (500.html)

**For each page, verify:**

- ✅ Page loads without errors
- ✅ JavaScript functionality works
- ✅ Styles are applied correctly
- ✅ No console errors (F12 → Console)
- ✅ HMR works (edit JS file, save, see instant update)

---

### Phase 4: Testing & Validation

#### Step 4.1: Update E2E Test Setup

**Edit `backend/tests_e2e/conftest.py`:**

```python
# Line 131: Update stats_file to manifest_file
manifest_file = Path(django.conf.settings.DJANGO_VITE["default"]["manifest_path"])

# Line 134: Update condition check
if _frontend_assets_need_rebuild(manifest_file):
    print("Building frontend assets for e2e tests...")
    subprocess.check_output(["pnpm", "--filter", "fpbase", "build"], stderr=subprocess.PIPE)

# Lines 138-142: Remove webpack-specific code
# No need to undo MockWebpackLoader (django-vite doesn't mock)

# Update _frontend_assets_need_rebuild function (lines 100-118)
def _frontend_assets_need_rebuild(manifest_file) -> bool:
    """Check if frontend assets need to be rebuilt."""
    if not manifest_file.is_file():
        return True

    # Check if manifest is valid JSON
    try:
        manifest = json.loads(manifest_file.read_bytes())
        if not manifest:
            return True
    except (json.JSONDecodeError, ValueError):
        return True

    # Check if source files are newer than manifest
    manifest_mtime = manifest_file.stat().st_mtime
    frontend_src = Path(__file__).parent.parent.parent / "frontend" / "src"
    if any(
        f.stat().st_mtime > manifest_mtime for f in frontend_src.rglob("*")
        if f.is_file() and not f.name.startswith(".")
    ):
        return True

    return False
```

#### Step 4.2: Run E2E Test Suite

```bash
cd backend

# Run e2e tests
uv run pytest backend/tests_e2e/ -v

# Expected: All tests pass
# If failures, investigate:
# - Are assets being built correctly?
# - Are pages loading correctly?
# - Are JavaScript interactions working?
```

#### Step 4.3: Run Unit Test Suite

```bash
cd backend

# Run unit tests
uv run pytest

# Expected: All tests pass
# Coverage should remain similar to baseline
```

---

### Phase 5: Production Deployment

#### Step 5.1: Remove Webpack Dependencies

**Edit `frontend/package.json`:**

```bash
cd frontend

# Remove webpack dependencies
pnpm remove webpack webpack-cli webpack-dev-server \
  webpack-bundle-tracker webpack-bundle-analyzer \
  clean-webpack-plugin copy-webpack-plugin \
  css-loader css-minimizer-webpack-plugin \
  mini-css-extract-plugin postcss-loader \
  sass-loader style-loader \
  babel-loader @babel/core @babel/preset-env @babel/preset-react \
  @babel/plugin-syntax-dynamic-import

# Remove webpack.config.js
rm webpack.config.js

# Update scripts (remove :old variants)
# Edit package.json:
{
  "scripts": {
    "dev": "vite",
    "build": "vite build",
    "preview": "vite preview"
  }
}
```

#### Step 5.2: Update CI/CD

**If using GitHub Actions, update `.github/workflows/*.yml`:**

```yaml
# Change build command
- name: Build frontend
  run: |
    cd frontend
    pnpm install
    pnpm run build  # Now uses Vite

# No changes needed for Django tests (they auto-detect manifest.json)
```

**If using Heroku, update `Procfile` (if needed):**

```
# Usually no changes needed
# Heroku buildpack should auto-detect pnpm and run build script
```

#### Step 5.3: Test Production Build

```bash
cd frontend

# Clean build
rm -rf dist

# Production build
NODE_ENV=production pnpm run build

# Verify output
ls -lah dist/

# Test with Django production settings
cd ../backend
DJANGO_DEBUG=False uv run python manage.py runserver

# Open http://localhost:8000
# Verify all pages load correctly
```

#### Step 5.4: Verify Sentry Source Maps

**After deploying to staging/production:**

1. Trigger an error (e.g., call undefined function in console)
2. Check Sentry for the error
3. Verify stack trace shows original source code (not minified)
4. Verify file paths are correct

#### Step 5.5: Update Documentation

**Edit `README.md`:**

```markdown
## Frontend Development

### Tech Stack

- **Build Tool:** Vite 7
- **Framework:** React 19
- **Languages:** JavaScript, Sass/SCSS

### Setup

```bash
cd frontend
pnpm install
```

### Development

```bash
# Start both Vite and Django dev servers concurrently
pnpm dev
```

Visit <http://localhost:8000>

### Production Build

```bash
pnpm run build
```

Output: `frontend/dist/`

### Entry Points

- `main`: Core site functionality (index.js)
- `embedscope`: Microscope embed viewer (embedscope.js)
- `blast`: BLAST search (blast-app.js)
- `spectraViewer`: Full spectra viewer (spectra-viewer.js)
- `simpleSpectraViewer`: Simple spectra widget (simple-spectra-viewer.js)
- `microscopeForm`: Microscope config form (microscope-form.js)
- `proteinTable`: Protein comparison table (protein-table.js)
- `litemol`: Protein structure viewer (my-litemol.js)

## Risk Assessment & Mitigation

### High Priority Risks

#### Risk 1: D3 v3/v7 Conflict in Embedscope

**Risk Level:** CRITICAL
**Probability:** Medium
**Impact:** High (chart breaks)

**Description:**
`embedscope.js` uses CDN D3 v3, but other bundles use D3 v7. If D3 v7 is accidentally bundled into embedscope, it will conflict with CDN D3 v3 and break charts.

**Mitigation:**

1. **Pre-migration:** Audit `src/embedscope.js` for D3 imports

   ```bash
   grep -r "from 'd3'" src/embedscope.js
   ```

   Expected: No matches

2. **Build-time check:** Add ESLint rule (see Phase 2.1)

3. **Manual verification:** After each build, inspect embedscope bundle:

   ```bash
   grep -i "d3" dist/assets/embedscope-*.js
   ```

   Expected: Only references to `window.d3` (CDN)

4. **Runtime test:** Load microscope embed page, check for console errors

**Rollback:** Revert vite.config.js changes, rebuild with webpack

---

#### Risk 2: Template Rendering Failures

**Risk Level:** HIGH
**Probability:** Medium
**Impact:** High (pages don't load)

**Description:**
10+ templates need migration from `{% render_bundle %}` to `{% vite_asset %}`. Incorrect paths or syntax will cause template errors.

**Mitigation:**

1. **Automated migration:** Use sed scripts (Phase 3.3)

2. **Manual verification:** Review each template change via git diff

3. **Comprehensive testing:** Test every page manually (Phase 3.5)

4. **E2E tests:** Run full test suite (Phase 4.2)

**Rollback:** Restore from `templates_backup_*.tar.gz`

---

#### Risk 3: Bundle Size Regression

**Risk Level:** MEDIUM
**Probability:** Low
**Impact:** Medium (slower page loads)

**Description:**
Vite's automatic code splitting may produce larger bundles than webpack's manual configuration.

**Mitigation:**

1. **Baseline measurement:** Record webpack bundle sizes before migration

2. **Vite profiling:** Compare Vite bundle sizes (Phase 2.2)

3. **Acceptance criteria:** ±10% variance acceptable, >20% requires investigation

4. **Optimization:** Adjust `manualChunks` in vite.config.js if needed

**Rollback:** Revert to webpack if bundles are >30% larger and optimization fails

---

### Medium Priority Risks

#### Risk 4: jQuery Global Injection

**Risk Level:** MEDIUM
**Probability:** Low
**Impact:** Medium (jQuery plugins break)

**Description:**
Webpack uses `ProvidePlugin` to auto-inject `$` and `jQuery`. Vite doesn't have this by default. Legacy code may expect global jQuery.

**Mitigation:**

1. **Manual injection:** In `src/index.js`:

   ```javascript
   import $ from 'jquery';
   window.$ = window.jQuery = $;
   ```

2. **Alternative:** Use `vite-plugin-inject` for auto-injection:

   ```javascript
   // vite.config.js
   import inject from '@rollup/plugin-inject';

   plugins: [
     inject({
       $: 'jquery',
       jQuery: 'jquery',
     }),
   ]
   ```

**Rollback:** None needed (fix is quick)

---

#### Risk 5: E2E Test Rebuild Detection

**Risk Level:** MEDIUM
**Probability:** Low
**Impact:** Medium (stale assets in tests)

**Description:**
E2E tests currently check `webpack-stats.json` freshness. Migration to `manifest.json` requires updating this logic.

**Mitigation:**

1. **Update `conftest.py`:** Change stats_file to manifest_file (Phase 4.1)

2. **Update `_frontend_assets_need_rebuild()`:** Check manifest.json validity

3. **Test:** Run e2e suite and verify rebuild triggers correctly

**Rollback:** Revert conftest.py changes

---

### Low Priority Risks

#### Risk 6: Sentry Source Map Upload

**Risk Level:** LOW
**Probability:** Low
**Impact:** Medium (degraded error tracking)

**Description:**
Sentry webpack plugin is replaced with Sentry Vite plugin. Configuration differences may break source map uploads.

**Mitigation:**

1. **Use official plugin:** `@sentry/vite-plugin`

2. **Test after deployment:** Trigger error, check Sentry for source-mapped stack

3. **Verify release:** Check Sentry releases for uploaded artifacts

**Rollback:** None needed (source maps are optional)

---

## Testing Strategy

### 1. Unit Tests

**No changes required** - Unit tests use MockWebpackLoader in test.py, which is replaced by django-vite's manifest-based loading.

**Verification:**

```bash
cd backend
uv run pytest --cov
```

**Success Criteria:**

- ✅ All tests pass
- ✅ Coverage remains ≥90% (or baseline)

---

### 2. E2E Tests

**Changes required** - Update `conftest.py` to use manifest.json instead of webpack-stats.json.

**Verification:**

```bash
cd backend
uv run pytest backend/tests_e2e/ -v
```

**Success Criteria:**

- ✅ All tests pass
- ✅ Frontend assets rebuild when source files change
- ✅ No stale asset issues

---

### 3. Manual Testing

**Comprehensive page testing** - Every page that uses JavaScript bundles.

**Test Matrix:**

| Page | Bundle(s) | Test Cases |
|------|----------|------------|
| Homepage | main | Navigation, search autocomplete, login |
| Protein detail | main, simpleSpectraViewer, litemol | Spectra widget, 3D viewer, favorites |
| Spectra viewer | spectraViewer | Load spectra, change filters, calculate efficiency |
| BLAST | blast | Submit search, view results |
| Protein table | proteinTable | Sort, filter, export |
| Microscope form | microscopeForm | Add/remove items, submit form |
| Microscope embed | main (css), embedscope | **D3 v3 chart renders correctly** |
| Compare | simpleSpectraViewer | Compare proteins, view spectra |
| Spectra graph widget | simpleSpectraViewer | Embedded spectra render |
| Error page (500) | main (css) | Styles load correctly |

**For each page:**

- ✅ Page loads without errors
- ✅ JavaScript functionality works
- ✅ Styles are applied correctly
- ✅ No console errors (F12 → Console)
- ✅ Network tab shows assets load from `/static/assets/` (production)
- ✅ Network tab shows assets load from `localhost:5173` (dev)

---

### 4. Performance Testing

**Metrics to measure:**

| Metric | Webpack (Baseline) | Vite (Target) | Change |
|--------|-------------------|---------------|--------|
| Production build time | 30-50s | 5-10s | -80% |
| Dev server startup | 10-30s | <1s | -90%+ |
| HMR update time | ~1000ms | ~100ms | -90% |
| Bundle size (main) | XXX KB | ±10% | Accept |
| Bundle size (embedscope) | XXX KB | ±10% | Accept |
| Bundle size (total) | XXX KB | ±10% | Accept |
| Lighthouse score | XX/100 | ≥XX/100 | Maintain |

**Tools:**

- Build time: `time pnpm run build`
- Bundle size: `du -sh dist/assets/*.js | sort -h`
- Lighthouse: Chrome DevTools → Lighthouse → Run audit

**Acceptance Criteria:**

- ✅ Build time improved by ≥50%
- ✅ HMR improved by ≥50%
- ✅ Bundle sizes within ±10% (or smaller)
- ✅ Lighthouse score unchanged or improved

---

## Reference Information

### Key Files Modified

**Frontend:**

- ✅ `frontend/vite.config.js` (new)
- ✅ `frontend/package.json` (scripts, dependencies)
- ✅ `frontend/src/index.js` (add modulepreload polyfill)
- ✅ `frontend/src/embedscope.js` (add modulepreload polyfill)
- ✅ `frontend/src/*.js` (all entry points: add modulepreload polyfill)
- ❌ `frontend/webpack.config.js` (deleted in Phase 5)

**Backend:**

- ✅ `backend/config/settings/base.py` (INSTALLED_APPS, DJANGO_VITE)
- ✅ `backend/config/settings/test.py` (DJANGO_VITE, remove MockWebpackLoader)
- ✅ `backend/tests_e2e/conftest.py` (_setup_frontend_assets,_frontend_assets_need_rebuild)

**Templates (10 files):**

- ✅ `backend/fpbase/templates/base.html`
- ✅ `backend/fpbase/templates/500.html`
- ✅ `backend/proteins/templates/spectra.html`
- ✅ `backend/proteins/templates/spectra_graph.html`
- ✅ `backend/proteins/templates/table.html`
- ✅ `backend/proteins/templates/compare.html`
- ✅ `backend/proteins/templates/proteins/microscope_form.html`
- ✅ `backend/proteins/templates/proteins/microscope_embed.html`
- ✅ `backend/proteins/templates/proteins/protein_detail.html`
- ✅ `backend/proteins/templates/proteins/blast.html`

---

### Documentation Resources

**Vite:**

- Official Guide: <https://vite.dev/guide/>
- Backend Integration: <https://vite.dev/guide/backend-integration.html>
- Migration from v6: <https://vite.dev/guide/migration>

**django-vite:**

- GitHub: <https://github.com/MrBin99/django-vite>
- PyPI: <https://pypi.org/project/django-vite/>
- README: <https://github.com/MrBin99/django-vite/blob/master/README.md>

**Sentry:**

- Vite Plugin: <https://docs.sentry.io/platforms/javascript/sourcemaps/uploading/vite/>

**Related Guides:**

- Webpack to Vite Migration: <https://www.sitepoint.com/webpack-vite-migration/>
- Vite vs Webpack 2025: <https://dualite.dev/blog/vite-vs-webpack>

---

### Commands Reference

**Vite:**

```bash
# Development
pnpm run dev                    # Start dev server (http://localhost:5173)
pnpm run build                  # Production build
pnpm run preview                # Preview production build

# Analysis
du -sh dist/assets/*.js         # Bundle sizes
```

**Django:**

```bash
# Development
uv run python manage.py runserver

# Testing
uv run pytest                   # All tests
uv run pytest backend/tests_e2e/ -v  # E2E tests only
```

**Combined Development:**

```bash
# Vite dev server and Django dev server concurrently
pnpm dev

# Visit: http://localhost:8000
```

---

### Migration Checklist

**Phase 1: Preparation**

- [ ] Node.js version ≥20.19 or ≥22.12
- [ ] Install Vite dependencies
- [ ] Create vite.config.js
- [ ] Add modulepreload polyfill to entry points
- [ ] Test standalone Vite build
- [ ] Verify manifest.json exists

**Phase 2: Configuration**

- [ ] Verify D3 v3/v7 isolation
- [ ] Profile bundle sizes (±10%)
- [ ] Test dev server + HMR
- [ ] Test Sentry integration

**Phase 3: Django Integration**

- [ ] Install django-vite
- [ ] Update settings (base.py, test.py)
- [ ] Migrate 10 templates
- [ ] Add HMR client to base.html
- [ ] Test each page manually

**Phase 4: Testing**

- [ ] Update E2E test setup (conftest.py)
- [ ] Run E2E test suite (all pass)
- [ ] Run unit test suite (all pass)
- [ ] Manual QA (all pages)

**Phase 5: Production**

- [ ] Remove webpack dependencies
- [ ] Update CI/CD scripts
- [ ] Test production build
- [ ] Verify Sentry source maps
- [ ] Update documentation
- [ ] Create PR and deploy
