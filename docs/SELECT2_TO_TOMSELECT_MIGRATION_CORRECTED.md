# Migration Plan: Select2 → Tom-Select & django-autocomplete-light → django-tomselect

**Status**: Ready for Implementation (Reviewed & Corrected by Architecture, Django, JavaScript, and Testing Experts)
**Created**: 2025-01-26
**Last Updated**: 2025-01-26
**Sentry Issue**: FPBASE-5H3 (JavaScript TypeError in Select2)
**Original Plan**: `SELECT2_TO_TOMSELECT_MIGRATION.md` (reference only)

---

## ⚠️ CRITICAL: Read This First

This is the **CORRECTED** migration guide incorporating expert feedback. Key corrections from original:

1. **jQuery CANNOT be removed** - Bootstrap 4 requires it (full removal needs Bootstrap 5 migration)
2. **Bundle savings: ~27% not 73%** - Realistic expectation: ~19KB per bundle
3. **django-tomselect API corrected** - Use `create_result_dict()` not non-existent methods
4. **Settings simplified** - No middleware/context processor needed
5. **Phase 0 added** - Pre-migration setup essential (baselines, testing infrastructure)
6. **embedscope.js** - Must keep Select2 temporarily (depends on legacy microscope.js)
7. **Deployment strategy added** - Feature flags, monitoring, rollout plan
8. **Testing enhanced** - TDD approach, accessibility automation, performance baselines

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Architectural Decision Record (ADR)](#architectural-decision-record-adr)
3. [Current State Audit](#current-state-audit)
4. [Target Architecture](#target-architecture)
5. [Migration Strategy](#migration-strategy)
6. [Detailed Migration Steps](#detailed-migration-steps)
   - [Phase 0: Pre-Migration Setup](#phase-0-pre-migration-setup)
   - [Phase 1: Simple Frontend Dropdowns](#phase-1-simple-frontend-dropdowns)
   - [Phase 2: Django Integration](#phase-2-django-integration)
   - [Phase 3: Custom AJAX](#phase-3-custom-ajax)
   - [Phase 4: Cleanup & Final Testing](#phase-4-cleanup--final-testing)
7. [Deployment Strategy](#deployment-strategy)
8. [Code Examples](#code-examples)
9. [Testing Strategy](#testing-strategy)
10. [Rollback Plan](#rollback-plan)
11. [Appendices](#appendices)
12. [Resources](#resources)

---

## Executive Summary

### Why Migrate?

**Current Issues**:

- Select2 4.0.13 is latest stable (from 2020, unmaintained)
- Select2 4.1.0-rc.0 has been RC for 5 years (no stable release)
- **Active Sentry issue FPBASE-5H3**: `TypeError: Cannot read properties of undefined (reading 'toString')`
- Larger bundle size (~30KB Select2 vs ~16KB Tom-Select)
- Missing modern accessibility improvements (ARIA, keyboard nav)
- Blocks Bootstrap 5 migration (technical debt)

**Benefits of Migration**:

- ✅ **Remove Select2 dependency** (~30KB saved)
- ✅ **~27% bundle size reduction** (~19KB net: -30KB Select2, +16KB Tom-Select, jQuery stays)
- ✅ **Fix Sentry issue FPBASE-5H3** (production error resolution)
- ✅ **Actively maintained** (Tom-Select last release: Feb 2025)
- ✅ **Better accessibility** (modern ARIA, keyboard navigation)
- ✅ **Future-proof** (enables Bootstrap 5 migration later)
- ✅ **Reduce custom SCSS** (remove most of 875-line Select2-Bootstrap theme)

**⚠️ IMPORTANT: jQuery Stays**

jQuery CANNOT be fully removed in this migration because:

- Bootstrap 4 JavaScript requires jQuery (modals, tooltips, popovers, etc.)
- `jquery.formset.js` used in microscope-form.js (dynamic formsets)
- 375+ jQuery calls throughout frontend/src/js/ codebase
- Full jQuery removal requires Bootstrap 5 migration (future ADR)

**Realistic Expectation**: This migration removes Select2 (~30KB), adds Tom-Select (~16KB), keeps jQuery (~30KB) = **~19KB net savings per bundle** (~27% reduction).

### Estimated Effort

| Phase | Effort | Risk Level | Blocking Issues |
|-------|--------|------------|-----------------|
| **Phase 0: Pre-Migration Setup** | 8-12 hours | Low | None - start immediately |
| **Phase 1: Simple Dropdowns** | 2-4 hours | Low | Requires Phase 0 complete |
| **Phase 2: Django Integration** | 1-2 days | Medium | Requires Phase 0, Phase 1 recommended |
| **Phase 3: Custom AJAX** | 1-2 hours | Low | Requires Phase 1 |
| **Phase 4: Cleanup & Testing** | 6-10 hours | Low | Requires all phases complete |
| **TOTAL** | **3-4 days** | **Low-Medium** | None if sequenced properly |

**Phase 0 includes**: Performance baselines, TDD test stubs, accessibility testing setup (axe-core), cross-browser config, feature flags, ADR documentation

---

## Architectural Decision Record (ADR)

### ADR-001: Migrate from Select2 to Tom-Select

**Date**: 2025-01-26
**Status**: Accepted
**Deciders**: Development Team (after expert reviews: Architecture, Django, JavaScript, Testing)

#### Context

FPbase uses Select2 4.0.13 (unmaintained since 2020) for autocomplete dropdowns. This creates:

1. **Production Errors**: Sentry FPBASE-5H3 - `TypeError` in Select2 code (`s.id.toString()` when `s.id` is undefined)
2. **Maintenance Risk**: No stable release in 5+ years, active development ceased
3. **Bundle Bloat**: 30KB Select2 + 30KB jQuery + 5KB theming = 65KB overhead
4. **Technical Debt**: Blocks Bootstrap 5 migration (needs jQuery removal)
5. **Accessibility Gaps**: Missing modern ARIA, keyboard nav improvements
6. **Developer Experience**: Poor docs, inactive community, no TypeScript support

**Alternatives Considered**:

- **Do Nothing**: Continue with Select2 → HIGH RISK (unmaintained, production errors)
- **Choices.js v11**: Actively maintained but larger bundle (41KB vs 16KB), less mature Django integration
- **Tom-Select**: Best choice - smaller, maintained, has django-tomselect package

#### Decision

Migrate to **Tom-Select 2.4.3** (frontend) + **django-tomselect 2025.9.1** (backend) using **incremental, phased approach** with **feature flags** for safe rollback.

**Migration Phases**:

1. **Phase 0**: Setup infrastructure (testing, baselines, feature flags)
2. **Phase 1**: Simple frontend dropdowns (microscope, FRET)
3. **Phase 2**: Django forms (lineage, spectrum, optical config, references)
4. **Phase 3**: Custom AJAX (project.js protein slug builder)
5. **Phase 4**: Cleanup (remove Select2, update tests, monitor production)

**Architectural Principles**:

- **Strangler Fig Pattern**: Old and new code coexist during migration
- **Incremental Migration**: Each phase independently testable/rollbackable
- **Feature Flag Safety**: `USE_TOMSELECT` env var for instant rollback
- **Test-Driven**: Write failing tests BEFORE code changes (red-green-refactor)

#### Consequences

**Positive**:

- ✅ Fixes FPBASE-5H3 production error
- ✅ ~19KB bundle reduction per affected bundle (27% savings)
- ✅ Future-proof (actively maintained, Feb 2025 release)
- ✅ Better accessibility (WCAG 2.1 AA compliant)
- ✅ Enables Bootstrap 5 migration (Tom-Select doesn't need jQuery)
- ✅ Better DX (modern API, good docs, TypeScript support)

**Negative**:

- ⚠️ 3-4 days engineering effort
- ⚠️ Team learning curve (new API, though well-documented)
- ⚠️ Migration risk (mitigated by testing, feature flags, rollback plan)
- ⚠️ jQuery stays (cannot remove until Bootstrap 5 migration)
- ⚠️ embedscope.js must keep Select2 (legacy microscope.js dependency)

**Mitigations**:

- Comprehensive testing (unit, integration, E2E, a11y, performance)
- Parallel implementation (old Select2 code stays until verified)
- Feature flag (`USE_TOMSELECT`) for instant production rollback
- TDD approach ensures no regressions
- Performance baselines verify improvements

#### Success Criteria

Migration succeeds when:

1. ✅ All tests passing (100% coverage maintained)
2. ✅ Bundle size reduced ≥15KB (target: ~19KB)
3. ✅ FPBASE-5H3 resolved (zero errors 7 days post-deploy)
4. ✅ WCAG 2.1 AA compliance maintained (axe-core clean scan)
5. ✅ Performance within ±5% baseline (page load, autocomplete latency)
6. ✅ Zero critical bugs first 2 weeks production
7. ✅ Successful rollback test in staging

#### Future Considerations

1. **Bootstrap 5 Migration** (6-12 months): Remove jQuery entirely
2. **React Component Wrapper** (future): Reusable Tom-Select React component
3. **embedscope.js Refactor** (future): Eliminate legacy microscope.js CDN dependency

---

## Current State Audit

### 1. Select2 Frontend Usage

#### Package Dependencies

**Location**: `frontend/package.json`

```json
{
  "dependencies": {
    "jquery": "^3.7.0",
    "select2": "^4.0.13",
    "select2-theme-bootstrap4": "^0.2.0-beta.6"
  },
  "devDependencies": {
    "bootstrap": "^4.6.2"
  }
}
```

#### Usage Patterns

**Pattern 1: Simple Dropdowns** (microscope-form.js:48-56, fret.js:113-122)

```javascript
$("#id_detector, #id_light_source, #id_extra_cameras, #id_extra_lights").select2({
  theme: "bootstrap",
  containerCssClass: ":all:",  // Copies CSS classes from <select>
  placeholder: "---------",
  allowClear: true,
  width: "auto",
})
```

**Pattern 2: AJAX Autocomplete** (project.js:207-228)

```javascript
$("#proteinSlug").select2({
  theme: "bootstrap",
  width: "80%",
  ajax: {
    url: "/autocomplete-protein",
    dataType: "json",
    cache: true,
    data: function(params) {
      return { q: params.term, type: "spectra", page: params.page }
    }
  }
})
```

**Pattern 3: Legacy CDN** (backend/fpbase/static/js/microscope.js:847)

- Loaded via CDN for backward compatibility
- Used on microscope detail pages
- **BLOCKER**: Cannot be migrated yet (see embedscope.js section)

**Pattern 4: Embedded Context** (frontend/src/embedscope.js:9)

```javascript
import "select2/dist/js/select2.full.js"
```

**⚠️ CRITICAL**: embedscope.js imports Select2 for legacy microscope.js (CDN-loaded). Must keep Select2 here until microscope.js refactored.

### 2. django-autocomplete-light Backend Usage

#### Package Dependencies

**Location**: `pyproject.toml`

```toml
dependencies = [
    "django-autocomplete-light>=3.11.0",
]
```

#### Autocomplete Views

**Location**: `backend/proteins/views/autocomplete.py`

All views extend `autocomplete.Select2QuerySetView`:

```python
from dal import autocomplete

class ProteinAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        qs = Protein.objects.with_spectra() if self.request.GET.get("type") == "spectra" else Protein.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs

    def get_results(self, context):
        return [{"id": r.slug, "text": r.name} for r in context["object_list"]]
```

**Similar views**: `LineageAutocomplete`, `StateAutocomplete`, `FilterAutocomplete`, `ReferenceAutocomplete`

#### Form Field Usage

**Files using DAL widgets**:

1. `backend/proteins/forms/forms.py:407-419` - LineageForm (parent field)
2. `backend/proteins/forms/spectrum.py:38-46` - SpectrumForm (owner_state field)
3. `backend/proteins/forms/microscope.py:309-324` - OpticalConfigForm (filter multi-select)
4. `backend/references/views.py:55-74` - ReferenceAutocomplete (requires auth)

### 3. Custom Styling

**Location**: `frontend/src/css/_select2-bootstrap.scss` (875 lines!)

Includes:

- Form validation states (`.is-invalid`, `.is-valid`)
- Input group integration
- Sizing classes (`-sm`, `-lg`)
- RTL support
- Dark theme support

**Strategy**: Remove most styling, but **keep validation state styles** for Django form compatibility.

### 4. End-to-End Tests

**Location**: `backend/fpbase/tests/test_end2end.py`

Tests verify:

- Select2 container initialization (`.select2-container`)
- Clickable state of widgets (`.select2-selection`)
- Results dropdown (`.select2-results__option`)
- Integration on microscope, FRET, spectra pages

**Selector Migration Required**:

| Select2 Selector | Tom-Select Selector | Usage |
|------------------|---------------------|-------|
| `.select2-container` | `.ts-wrapper` | Container element |
| `.select2-selection` | `.ts-control` | Click target |
| `#select2-{id}-container` | `#{id} + .ts-wrapper` | Specific instance |
| `.select2-results` | `.ts-dropdown` | Dropdown container |
| `.select2-results__option` | `.ts-dropdown .option` | Individual options |

---

## Target Architecture

### Tom-Select (Frontend)

**Package**: `tom-select` v2.4.3 (Feb 2025 release)

**Features**:

- Framework-agnostic (no jQuery dependency)
- ~16KB gzipped (vs Select2's ~30KB)
- Smart search with diacritics support
- AJAX/remote data loading
- Multi-select, tagging, item creation
- Modern accessibility (ARIA, keyboard nav)
- Plugin system (remove_button, clear_button, dropdown_header, etc.)

**Basic Initialization**:

```javascript
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

new TomSelect('#select-element', {
  placeholder: 'Select an option...',
  maxItems: 1,  // Single select (null = multi-select)
  plugins: ['clear_button'],

  // AJAX
  load: function(query, callback) {
    fetch(`/api/search?q=${encodeURIComponent(query)}`)
      .then(response => response.json())
      .then(json => callback(json.results))
      .catch(() => callback())
  },

  // Configuration
  valueField: 'id',
  labelField: 'text',
  searchField: 'text',
})
```

### django-tomselect (Backend)

**Package**: `django-tomselect` v2025.9.1 (Sep 2025 release)

**Features**:

- Drop-in replacement for django-autocomplete-light
- Built on Tom-Select (no jQuery)
- Bootstrap 4/5 theming
- Permission handling, security, i18n
- Dependent field filtering (`filter_by`)
- Plugin support (clear_button, remove_button, etc.)

**View Pattern** (✅ CORRECTED):

```python
from django_tomselect.autocompletes import AutocompleteModelView

class ProteinTomSelectView(AutocompleteModelView):
    model = Protein
    search_lookups = ["name__icontains"]
    ordering = ["name"]
    page_size = 20

    # ✅ CORRECT: Use create_result_dict() to customize id/text
    def create_result_dict(self, result):
        """Customize result dictionary (replaces get_result_value/label)."""
        return {
            "id": result.slug,  # Use slug instead of pk
            "text": result.name,
        }

    def get_queryset(self):
        """Apply custom filtering."""
        qs = super().get_queryset()
        if self.request.GET.get("type") == "spectra":
            qs = qs.filter(default_state__spectra__isnull=False).distinct()
        return qs
```

**Form Pattern**:

```python
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginClearButton

class MyForm(forms.ModelForm):
    protein = TomSelectModelChoiceField(
        config=TomSelectConfig(
            url="protein-autocomplete",
            placeholder="Select a protein...",
            minimum_query_length=2,
            highlight=True,
            preload="focus",
            plugin_clear_button=PluginClearButton(title="Clear"),
            css_framework="bootstrap4",
        ),
        queryset=Protein.objects.all(),
        required=False,
    )
```

---

## Migration Strategy

### Approach: Incremental, Test-Driven, Low-Risk

**Principles**:

1. ✅ **TDD**: Write failing tests BEFORE code changes (red-green-refactor)
2. ✅ **Incremental**: Migrate in phases, not all-at-once
3. ✅ **Parallel Implementation**: Keep old code until new verified
4. ✅ **Feature Flags**: `USE_TOMSELECT` env var for instant rollback
5. ✅ **Measure Everything**: Baselines for bundle size, performance, errors

### Phase Overview

#### Phase 0: Pre-Migration Setup (8-12 hours) - NEW

**Goal**: Establish infrastructure before touching production code

**Tasks**:
1. Measure performance baselines (bundle size, page load, autocomplete latency)
2. Set up accessibility testing (install axe-core, create test suite)
3. Configure cross-browser testing (Playwright multi-browser)
4. Create TDD test stubs (write failing tests for Tom-Select)
5. Implement feature flag system (`USE_TOMSELECT`)
6. Set up Sentry monitoring (alerts, error budgets)
7. Document ADR (this document)

**Deliverables**:
- `baseline_metrics.json` (bundle sizes, performance data)
- `backend/fpbase/tests/test_tomselect_migration.py` (TDD stubs marked `@pytest.mark.xfail`)
- `backend/fpbase/tests/test_accessibility.py` (axe-core integration)
- Feature flag in settings: `USE_TOMSELECT = env.bool("USE_TOMSELECT", default=False)`
- Sentry alerts configured
- ADR documented

#### Phase 1: Simple Frontend Dropdowns (2-4 hours)

**Goal**: Replace simple Select2 dropdowns with Tom-Select

**Scope**:
- `microscope-form.js` - camera/light dropdowns
- `fret.js` - donor/acceptor selection

**Why First**: Lowest risk, fastest feedback, no backend changes, builds team confidence

**TDD Approach**:
1. Run failing test: `test_microscope_form_uses_tomselect()` (should fail)
2. Migrate code
3. Run test: should pass
4. Refactor for quality
5. Run full suite: all pass

#### Phase 2: Django Integration (1-2 days)

**Goal**: Replace django-autocomplete-light with django-tomselect

**Scope**:
- Install django-tomselect
- Create new autocomplete views (parallel to DAL)
- Migrate forms one-by-one (simplest → complex)
- Update URL patterns
- Test thoroughly

**Migration Order**:
1. LineageForm (single select, simple)
2. SpectrumForm (single select with custom display)
3. OpticalConfigForm (multi-select)
4. ReferenceAutocomplete (requires authentication)

#### Phase 3: Custom AJAX (1-2 hours)

**Goal**: Replace Select2 AJAX in project.js

**Scope**: Protein slug builder autocomplete

#### Phase 4: Cleanup & Final Testing (6-10 hours)

**Goal**: Remove old dependencies, comprehensive testing

**Tasks**:
1. Remove Select2/theme from package.json
2. Remove django-autocomplete-light from pyproject.toml
3. Remove most of _select2-bootstrap.scss (keep validation styles)
4. Remove old DAL views/URLs
5. Update E2E tests (all Select2 selectors → Tom-Select)
6. Run full test suite
7. Performance verification
8. Accessibility scan
9. Cross-browser testing
10. Monitor Sentry for errors

---

## Detailed Migration Steps

### Phase 0: Pre-Migration Setup

#### Step 0.1: Measure Performance Baselines

Create baseline collection script:

```bash
# scripts/collect_baselines.sh

#!/bin/bash
set -e

echo "📊 Collecting Performance Baselines..."

# Bundle sizes
echo "Bundle sizes (before migration):" > baseline_metrics.json
pnpm build
du -b backend/fpbase/static/bundles/*.js | awk '{sum+=$1} END {print "Total: " sum " bytes"}' >> baseline_metrics.json

# Page load times (requires running server)
echo "Starting dev server..."
uv run backend/manage.py runserver &
SERVER_PID=$!
sleep 5

# Use Lighthouse CLI or curl timing
curl -w "@curl-format.txt" -o /dev/null -s "http://localhost:8000/microscopes/add/" >> baseline_metrics.json

kill $SERVER_PID

echo "✅ Baselines saved to baseline_metrics.json"
cat baseline_metrics.json
```

Run it:

```bash
chmod +x scripts/collect_baselines.sh
./scripts/collect_baselines.sh
```

#### Step 0.2: Set Up Accessibility Testing

Install axe-core:

```bash
pnpm add --save-dev @axe-core/playwright axe-core
```

Create accessibility test file:

```python
# backend/fpbase/tests/test_accessibility.py

import pytest
from axe_playwright_python import Axe

@pytest.fixture
def axe(page):
    """Provide axe-core accessibility scanner."""
    return Axe(page)

@pytest.mark.accessibility
class TestAccessibility:
    """WCAG 2.1 AA compliance tests."""

    def test_microscope_form_wcag_compliance(self, page, live_server, axe):
        """Microscope form must meet WCAG 2.1 AA standards."""
        page.goto(f'{live_server.url}/microscopes/add/')

        # Run accessibility scan
        results = axe.run()

        # Assert no violations
        violations = results['violations']
        assert len(violations) == 0, \
            f"WCAG violations found: {axe.report(violations)}"
```

Run baseline scan:

```bash
uv run pytest backend/fpbase/tests/test_accessibility.py -v
# Should pass (or document existing violations to fix separately)
```

#### Step 0.3: Create TDD Test Stubs

Create migration test file:

```python
# backend/fpbase/tests/test_tomselect_migration.py

import pytest

class TestTomSelectMigration:
    """TDD tests for Select2 → Tom-Select migration.

    These tests are marked xfail initially and should PASS after migration.
    """

    @pytest.mark.xfail(reason="Tom-Select not yet implemented")
    def test_microscope_form_uses_tomselect(self, page, live_server):
        """Microscope form should use Tom-Select, not Select2."""
        page.goto(f'{live_server.url}/microscopes/add/')

        # Should find Tom-Select wrapper
        ts_wrapper = page.query_selector('.ts-wrapper')
        assert ts_wrapper is not None, "Tom-Select not initialized"

        # Should NOT find Select2 container
        select2 = page.query_selector('.select2-container')
        assert select2 is None, "Select2 still present"

    @pytest.mark.xfail(reason="django-tomselect not yet implemented")
    def test_lineage_form_uses_tomselect_field(self):
        """LineageForm.parent should use TomSelectModelChoiceField."""
        from proteins.forms.forms import LineageForm
        from django_tomselect.forms import TomSelectModelChoiceField

        form = LineageForm()
        assert isinstance(form.fields['parent'], TomSelectModelChoiceField)
```

Run to verify they fail:

```bash
uv run pytest backend/fpbase/tests/test_tomselect_migration.py -v
# Should show 2 xfail (expected failures)
```

#### Step 0.4: Implement Feature Flag

Add to settings:

```python
# backend/config/settings/base.py

# Feature flag for Tom-Select migration (rollback safety)
USE_TOMSELECT = env.bool("USE_TOMSELECT", default=False)
```

Use in code (example for Phase 2):

```python
# backend/proteins/forms/forms.py

from django.conf import settings

if settings.USE_TOMSELECT:
    from django_tomselect.forms import TomSelectModelChoiceField
    AutocompleteField = TomSelectModelChoiceField
else:
    from dal.autocomplete import ModelSelect2
    # Wrapper for compatibility
    class AutocompleteField:
        pass  # Keep DAL implementation
```

#### Step 0.5: Configure Sentry Monitoring

```python
# backend/config/settings/production.py

import sentry_sdk

# Configure release-specific alerts
SENTRY_ENVIRONMENT = env("SENTRY_ENVIRONMENT", default="production")

sentry_sdk.init(
    # ... existing config ...

    # Tag events for migration monitoring
    before_send=lambda event, hint: {
        **event,
        "tags": {
            **event.get("tags", {}),
            "tomselect_migration": "active" if USE_TOMSELECT else "inactive",
        }
    }
)
```

Set up Sentry alert:

1. Go to Sentry → Alerts → Create Alert
2. Filter: `tags.tomselect_migration:active AND message contains "TomSelect"`
3. Threshold: >10 events in 1 hour
4. Action: Email team, post to #engineering-alerts

#### Step 0.6: Configure Cross-Browser Testing

Update Playwright config:

```javascript
// playwright.config.js

module.exports = {
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
    },
    {
      name: 'firefox',
      use: { ...devices['Desktop Firefox'] },
    },
    {
      name: 'webkit',
      use: { ...devices['Desktop Safari'] },
    },
  ],
}
```

#### Step 0.7: Verification Checklist

Before proceeding to Phase 1:

- [ ] Baseline metrics collected (`baseline_metrics.json` exists)
- [ ] Accessibility tests passing (or violations documented)
- [ ] TDD stubs created and failing (xfail status)
- [ ] Feature flag implemented (`USE_TOMSELECT` in settings)
- [ ] Sentry alerts configured
- [ ] Cross-browser config ready
- [ ] Team briefed on migration plan

---

### Phase 1: Simple Frontend Dropdowns

#### Step 1.1: Install Tom-Select

```bash
pnpm add tom-select
```

#### Step 1.2: Update microscope-form.js (✅ CORRECTED)

**Before**:

```javascript
import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
import "select2/dist/js/select2.full"

$("#id_detector, #id_light_source, #id_extra_cameras, #id_extra_lights").select2({
  theme: "bootstrap",
  containerCssClass: ":all:",
  placeholder: "---------",
  allowClear: true,
  width: "auto",
})
```

**After** (✅ handles containerCssClass and width):

```javascript
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

const selectElements = [
  '#id_detector',
  '#id_light_source',
  '#id_extra_cameras',
  '#id_extra_lights'
]

selectElements.forEach(selector => {
  const element = document.querySelector(selector)
  if (!element) return

  new TomSelect(selector, {
    plugins: ['clear_button'],
    placeholder: '---------',
    allowEmptyOption: true,

    // ✅ Copy CSS classes from original select (matches containerCssClass: ":all:")
    controlClass: element.className,

    // ✅ Tom-Select defaults to 100% width; for auto behavior, set explicitly
    // (Tom-Select doesn't have width: "auto" - it's always full-width of container)
  })
})
```

**Note**: `width: "auto"` behavior differs - Tom-Select is always 100% of container. Ensure container has correct width.

#### Step 1.3: Update fret.js

**Before**:

```javascript
$("#donor-select, #acceptor-select").select2({
  theme: "bootstrap",
  containerCssClass: ":all:",
  width: "auto",
})
```

**After**:

```javascript
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

['#donor-select', '#acceptor-select'].forEach(selector => {
  const element = document.querySelector(selector)
  if (!element) return

  new TomSelect(selector, {
    controlClass: element.className,  // Match containerCssClass: ":all:"
  })
})
```

#### Step 1.4: Run TDD Tests

```bash
# This test should now PASS (remove xfail marker first)
uv run pytest backend/fpbase/tests/test_tomselect_migration.py::TestTomSelectMigration::test_microscope_form_uses_tomselect -v

# Should output: PASSED ✅
```

#### Step 1.5: Manual Testing

Start dev server:

```bash
pnpm dev
```

Test in browser:

1. Navigate to `/microscopes/add/`
2. Click detector dropdown - should open Tom-Select dropdown
3. Select an option - should populate
4. Click clear button (X) - should clear selection
5. Navigate to `/fret/` calculator
6. Test donor/acceptor selection
7. Verify styling matches Bootstrap 4 theme

#### Step 1.6: Verify Bundle Size Reduction

```bash
pnpm build

# Check bundle size
du -h backend/fpbase/static/bundles/microscope-form.*.js
du -h backend/fpbase/static/bundles/fret.*.js

# Compare against baseline_metrics.json
# Should see ~19KB reduction per bundle
```

#### Step 1.7: Update embedscope.js (⚠️ SPECIAL HANDLING)

**CRITICAL**: embedscope.js CANNOT migrate yet because legacy microscope.js (CDN-loaded) depends on Select2.

**Strategy**: Keep Select2 in embedscope.js, add Tom-Select for new features:

```javascript
// frontend/src/embedscope.js

// ✅ KEEP Select2 for legacy microscope.js compatibility
import "select2/dist/js/select2.full.js"

// ✅ ADD Tom-Select for new features
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

// Use Tom-Select for new widgets, Select2 for legacy
```

**Document as tech debt**:

```markdown
## TODO: embedscope.js Select2 Removal

**Blocker**: backend/fpbase/static/js/microscope.js (847 lines) loaded via CDN expects Select2

**Required Before Removal**:
1. Refactor microscope.js to use Tom-Select OR
2. Bundle microscope.js with webpack (no CDN) OR
3. Deprecate microscope.js entirely (rewrite functionality)

**Estimated Effort**: 2-3 days
**Priority**: Low (not blocking migration)
**Tracked In**: GitHub Issue #XXX
```

---

### Phase 2: Django Integration

#### Step 2.1: Install django-tomselect

```bash
uv add django-tomselect
```

#### Step 2.2: Configure Django Settings (✅ CORRECTED)

**File**: `backend/config/settings/base.py`

```python
# ✅ CORRECT: Only add to INSTALLED_APPS (no middleware/context processor needed!)
INSTALLED_APPS = [
    # ... existing apps
    "dal",  # Keep during migration for parallel implementation
    "dal_select2",
    "django_tomselect",  # Add new
]

# Optional: Global Tom-Select configuration
TOMSELECT_CONFIG = {
    "css_framework": "bootstrap4",  # Match current Bootstrap 4
}
```

**❌ WRONG** (from original doc - ignore these):

```python
# ❌ DO NOT ADD - These don't exist in django-tomselect!
MIDDLEWARE = [
    "django_tomselect.middleware.TomSelectMiddleware",  # ❌ Doesn't exist
]

TEMPLATES = [{
    "context_processors": [
        "django_tomselect.context_processors.tomselect",  # ❌ Doesn't exist
    ],
}]
```

#### Step 2.3: Collect Static Files

```bash
uv run backend/manage.py collectstatic --noinput
```

#### Step 2.4: Create Autocomplete Views (✅ CORRECTED)

**File**: `backend/proteins/views/autocomplete.py`

Add new Tom-Select views alongside existing DAL views:

```python
from __future__ import annotations

from django.contrib.auth.mixins import LoginRequiredMixin
from django_tomselect.autocompletes import AutocompleteModelView

from ..models import Filter, Lineage, Protein, State


class ProteinTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete for Protein model."""

    model = Protein
    search_lookups = ["name__icontains", "slug__icontains"]
    ordering = ["name"]
    page_size = 20

    # ✅ CORRECT: Use create_result_dict() to customize response
    def create_result_dict(self, result):
        """Return custom result with slug as ID (not pk)."""
        return {
            "id": result.slug,
            "text": result.name,
        }

    def get_queryset(self):
        """Apply type-based filtering."""
        qs = super().get_queryset()
        if self.request.GET.get("type") == "spectra":
            qs = qs.filter(default_state__spectra__isnull=False).distinct()
        return qs


class LineageTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete for Lineage model."""

    model = Lineage
    search_lookups = ["protein__name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        """Display protein name for lineage selection."""
        return {
            "id": result.pk,
            "text": result.protein.name,  # ✅ Access relationship in Python, not config
        }


class StateTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete for State model."""

    model = State
    search_lookups = ["protein__name__icontains", "name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        """Use State's __str__ method for display."""
        return {
            "id": result.pk,
            "text": str(result),  # ✅ Calls State.__str__()
        }


class FilterTomSelectView(AutocompleteModelView):
    """Tom-Select autocomplete for Filter model."""

    model = Filter
    search_lookups = ["name__icontains", "part__icontains"]
    ordering = ["name"]
    page_size = 20
    # No create_result_dict needed - defaults to pk/str(obj) work fine


# ✅ KEEP old DAL views during migration (parallel implementation)
from dal import autocomplete  # noqa: E402


class ProteinAutocomplete(autocomplete.Select2QuerySetView):
    """DEPRECATED: Use ProteinTomSelectView instead. Remove after migration."""

    def get_results(self, context):
        return [
            {"id": result.slug, "text": result.name}
            for result in context["object_list"]
        ]

    def get_queryset(self):
        qs = Protein.objects.with_spectra() if self.request.GET.get("type") == "spectra" else Protein.objects.all()
        if self.q:
            qs = qs.filter(name__icontains=self.q)
        return qs


# ... similar for LineageAutocomplete, StateAutocomplete, FilterAutocomplete
```

#### Step 2.5: Add URL Patterns

**File**: `backend/proteins/urls.py`

```python
from .views.autocomplete import (
    # Old DAL views (keep for backward compat during migration)
    ProteinAutocomplete,
    LineageAutocomplete,
    StateAutocomplete,
    FilterAutocomplete,
    # New Tom-Select views
    ProteinTomSelectView,
    LineageTomSelectView,
    StateTomSelectView,
    FilterTomSelectView,
)

urlpatterns = [
    # ... existing patterns

    # Old DAL endpoints (DEPRECATED - remove in Phase 4)
    path("autocomplete-protein/", ProteinAutocomplete.as_view(), name="protein-autocomplete"),
    path("lineage-autocomplete/", LineageAutocomplete.as_view(), name="lineage-autocomplete"),
    path("state-autocomplete/", StateAutocomplete.as_view(), name="state-autocomplete"),
    path("filter-autocomplete/", FilterAutocomplete.as_view(), name="filter-autocomplete"),

    # New Tom-Select endpoints (use 'ts/' prefix for clarity)
    path("ts/protein/", ProteinTomSelectView.as_view(), name="protein-ts"),
    path("ts/lineage/", LineageTomSelectView.as_view(), name="lineage-ts"),
    path("ts/state/", StateTomSelectView.as_view(), name="state-ts"),
    path("ts/filter/", FilterTomSelectView.as_view(), name="filter-ts"),
]
```

#### Step 2.6: Migrate LineageForm (Simplest First)

**File**: `backend/proteins/forms/forms.py`

**Before**:

```python
from dal import autocomplete

class LineageForm(forms.ModelForm):
    parent = forms.ModelChoiceField(
        queryset=Lineage.objects.all().prefetch_related("protein"),
        widget=autocomplete.ModelSelect2(
            url="proteins:lineage-autocomplete",
            attrs={
                "data-theme": "bootstrap",
                "data-width": "100%",
                "data-placeholder": "----------",
            },
        ),
        required=False,
        help_text="Direct evolutionary ancestor",
    )
```

**After**:

```python
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginClearButton

class LineageForm(forms.ModelForm):
    parent = TomSelectModelChoiceField(
        config=TomSelectConfig(
            url="proteins:lineage-ts",  # ✅ Use new Tom-Select URL
            placeholder="Select parent lineage",
            minimum_query_length=2,
            highlight=True,
            preload="focus",
            plugin_clear_button=PluginClearButton(title="Clear Selection"),
            css_framework="bootstrap4",
        ),
        queryset=Lineage.objects.select_related("protein"),
        required=False,
        help_text="Direct evolutionary ancestor (must have ancestor of its own to appear)",
    )

    class Meta:
        model = Lineage
        fields = ("protein", "parent", "mutation")

    # ... existing clean() and __init__() methods stay the same
```

#### Step 2.7: Test LineageForm

**TDD test should now pass**:

```bash
uv run pytest backend/fpbase/tests/test_tomselect_migration.py::test_lineage_form_uses_tomselect_field -v
# Should PASS (remove xfail marker first)
```

**Integration test**:

```bash
uv run pytest backend/proteins/tests/test_forms.py::TestLineageForm -v
```

**E2E test**:

```python
# backend/fpbase/tests/test_end2end.py

def test_lineage_form_autocomplete_e2e(page, live_server):
    """Test lineage parent selection with Tom-Select."""
    page.goto(f'{live_server.url}/admin/proteins/lineage/add/')

    # Wait for Tom-Select initialization
    parent_wrapper = page.wait_for_selector('.ts-wrapper', timeout=5000)

    # Type to search
    parent_input = parent_wrapper.locator('input[type="text"]')
    parent_input.type('GFP', delay=100)

    # Wait for results
    page.wait_for_load_state('networkidle')
    first_option = page.wait_for_selector('.ts-dropdown .option:first-child', state='visible')
    first_option.click()

    # Verify selection
    selected_item = parent_wrapper.locator('.item')
    expect(selected_item).to_contain_text('GFP')
```

#### Step 2.8: Migrate SpectrumForm

**File**: `backend/proteins/forms/spectrum.py`

```python
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig

class SpectrumForm(forms.ModelForm):
    owner_state = TomSelectModelChoiceField(
        config=TomSelectConfig(
            url="proteins:state-ts",
            placeholder="Select a state...",
            minimum_query_length=2,
            highlight=True,
            preload="focus",
            css_framework="bootstrap4",
        ),
        queryset=State.objects.select_related("protein"),
        required=False,
        label=mark_safe('Protein<span class="asteriskField">*</span>'),
    )
```

#### Step 2.9: Migrate OpticalConfigForm (Multi-Select)

**File**: `backend/proteins/forms/microscope.py`

```python
from django_tomselect.forms import TomSelectModelMultipleChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginRemoveButton

class MultipleFilterField(TomSelectModelMultipleChoiceField):
    def __init__(self, label):
        super().__init__(
            label=label,
            config=TomSelectConfig(
                url="proteins:filter-ts",
                placeholder="----------",
                minimum_query_length=1,
                highlight=True,
                preload="focus",
                max_items=None,  # ✅ Unlimited multi-select
                plugin_remove_button=PluginRemoveButton(title="Remove"),
                css_framework="bootstrap4",
            ),
            queryset=Filter.objects.all(),
            required=False,
        )

class OpticalConfigForm(forms.ModelForm):
    ex_filters = MultipleFilterField("Excitation Filter(s)")
    em_filters = MultipleFilterField("Emission Filter(s)")
    ref_em_filters = MultipleFilterField("Advanced: Reflective Emission Filter(s)")
    bs_filters = MultipleFilterField("Dichroic Filter(s)")
```

#### Step 2.10: Migrate ReferenceAutocomplete (With Auth)

**File**: `backend/references/views.py`

```python
from django.contrib.auth.mixins import LoginRequiredMixin
from django_tomselect.autocompletes import AutocompleteModelView

class ReferenceTomSelectView(LoginRequiredMixin, AutocompleteModelView):
    """Autocomplete for references - requires authentication."""

    model = Reference
    search_lookups = ["doi__icontains"]
    ordering = ["citation"]
    page_size = 20

    # ✅ Use standard Django LoginRequiredMixin (not permission_required attribute)

    def create_result_dict(self, result):
        return {
            "id": result.doi,
            "text": result.citation,
        }
```

#### Step 2.11: Verify Template Includes form.media

Check all form templates include media:

```bash
grep -r "{{ form.media }}" backend/fpbase/templates/
# Should find usage in form templates
```

If missing, add to templates:

```django
{# In form template #}
{% block extra_css %}
    {{ form.media.css }}
{% endblock %}

{% block extra_js %}
    {{ form.media.js }}
{% endblock %}
```

---

### Phase 3: Custom AJAX

#### Step 3.1: Update project.js

**File**: `frontend/src/js/project.js`

**Before**:

```javascript
$("#proteinSlug").select2({
  theme: "bootstrap",
  width: "80%",
  ajax: {
    url: "/autocomplete-protein",
    dataType: "json",
    cache: true,
    data: function(params) {
      return {
        q: params.term,
        type: "spectra",
        page: params.page,
      }
    }
  }
})
```

**After** (✅ with error handling):

```javascript
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

new TomSelect('#proteinSlug', {
  valueField: 'id',
  labelField: 'text',
  searchField: 'text',

  load: function(query, callback) {
    if (!query.length) return callback()

    const url = `/proteins/ts/protein/?q=${encodeURIComponent(query)}&type=spectra`

    fetch(url)
      .then(response => {
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}: ${response.statusText}`)
        }
        return response.json()
      })
      .then(json => {
        callback(json.results || [])
      })
      .catch(error => {
        console.error('Autocomplete fetch failed:', error)
        // ✅ Show user feedback
        this.settings.placeholder = 'Search temporarily unavailable'
        callback()  // Empty results on error
      })
  },

  render: {
    option: function(item, escape) {
      return `<div>${escape(item.text)}</div>`
    }
  }
})
```

#### Step 3.2: Test

```bash
pnpm dev

# Test in browser:
# 1. Navigate to protein slug builder page
# 2. Type in autocomplete
# 3. Verify results appear
# 4. Select a protein
# 5. Verify URL building still works
```

---

### Phase 4: Cleanup & Final Testing

#### Step 4.1: Remove Old Frontend Dependencies

```bash
pnpm remove select2 select2-theme-bootstrap4
```

**✅ Keep jQuery** (still needed for Bootstrap 4 + jquery.formset):

```bash
# DO NOT remove jQuery!
# pnpm remove jquery  # ❌ DON'T DO THIS
```

#### Step 4.2: Remove Old Backend Dependencies

```bash
uv remove django-autocomplete-light
uv sync
```

#### Step 4.3: Remove Old Code

**Delete old DAL views**:

```python
# In backend/proteins/views/autocomplete.py
# Delete these classes (keep only Tom-Select versions):
# - ProteinAutocomplete
# - LineageAutocomplete
# - StateAutocomplete
# - FilterAutocomplete
```

**Delete old URL patterns**:

```python
# In backend/proteins/urls.py
# Remove these paths:
# path("autocomplete-protein/", ...)
# path("lineage-autocomplete/", ...)
# path("state-autocomplete/", ...)
# path("filter-autocomplete/", ...)
```

**Remove Select2 imports**:

```bash
# Find and remove
grep -r "import.*select2" frontend/src/ --include="*.js"
# Except in embedscope.js (must keep for legacy microscope.js)
```

#### Step 4.4: Update SCSS (Partial Removal)

**File**: `frontend/src/css/_select2-bootstrap.scss`

**Strategy**: Remove most styling, but **keep validation state styles** for Django forms:

```scss
// ✅ KEEP: Tom-Select validation state styles for Django
.ts-wrapper.is-invalid {
  border-color: #dc3545;

  .ts-control {
    border-color: #dc3545;
  }
}

.ts-wrapper.is-valid {
  border-color: #28a745;

  .ts-control {
    border-color: #28a745;
  }
}

// ❌ DELETE: All other Select2-Bootstrap styles (input groups, sizing, etc.)
```

Rename file:

```bash
mv frontend/src/css/_select2-bootstrap.scss frontend/src/css/_tomselect-django.scss
```

Update imports in main SCSS file.

#### Step 4.5: Update E2E Tests (All Select2 Selectors)

**File**: `backend/fpbase/tests/test_end2end.py`

Find all Select2 selectors and replace:

```python
# Before
select2_container = page.query_selector('.select2-container')

# After
ts_wrapper = page.query_selector('.ts-wrapper')

# Before
page.click('.select2-selection')

# After
page.click('.ts-control')

# Before
results = page.query_selector_all('.select2-results__option')

# After
results = page.query_selector_all('.ts-dropdown .option')
```

**Complete selector migration table** (for reference):

| Old (Select2) | New (Tom-Select) | Usage |
|---------------|------------------|-------|
| `.select2-container` | `.ts-wrapper` | Container element |
| `.select2-selection` | `.ts-control` | Click target |
| `.select2-search__field` | `.ts-control input[type="text"]` | Search input |
| `#select2-{id}-container` | `#{id} + .ts-wrapper` | Specific instance |
| `.select2-results` | `.ts-dropdown` | Results container |
| `.select2-results__option` | `.ts-dropdown .option` | Individual options |
| `.select2-selection__clear` | `.clear-button` | Clear button |
| `.select2-selection__choice__remove` | `.remove` | Remove tag (multi) |

#### Step 4.6: Run Full Test Suite

```bash
# Python tests
uv run pytest --cov --cov-report=html --cov-fail-under=95

# Should output: >95% coverage, all tests passing
```

#### Step 4.7: Verify Bundle Sizes

```bash
pnpm build

# Compare against baseline
du -h backend/fpbase/static/bundles/*.js | tee bundle_sizes_after.txt

# Compare
diff baseline_metrics.json bundle_sizes_after.txt
# Should show ~19KB reduction per bundle
```

#### Step 4.8: Run Accessibility Scan

```bash
uv run pytest backend/fpbase/tests/test_accessibility.py -v

# Should PASS with zero WCAG violations
```

#### Step 4.9: Cross-Browser Testing

```bash
# Run E2E tests across all browsers
uv run pytest backend/fpbase/tests/test_end2end.py --browser chromium -v
uv run pytest backend/fpbase/tests/test_end2end.py --browser firefox -v
uv run pytest backend/fpbase/tests/test_end2end.py --browser webkit -v

# All should PASS
```

#### Step 4.10: Performance Testing

Create performance test script:

```python
# backend/fpbase/tests/test_performance.py

import pytest
import time

@pytest.mark.performance
class TestTomSelectPerformance:

    def test_dropdown_open_speed(self, page, live_server):
        """Dropdown should open in <50ms."""
        page.goto(f'{live_server.url}/microscopes/add/')

        start = time.perf_counter()
        page.click('.ts-control')
        page.wait_for_selector('.ts-dropdown[style*="display"]', timeout=1000)
        duration = (time.perf_counter() - start) * 1000

        assert duration < 50, f"Dropdown took {duration}ms (limit: 50ms)"

    def test_autocomplete_latency(self, page, live_server):
        """Autocomplete should respond in <200ms."""
        page.goto(f'{live_server.url}/proteins/lineage/add/')

        input_field = page.locator('.ts-wrapper input')

        start = time.perf_counter()
        input_field.type('GFP')
        page.wait_for_selector('.ts-dropdown .option', timeout=2000)
        duration = (time.perf_counter() - start) * 1000

        assert duration < 200, f"Autocomplete took {duration}ms (limit: 200ms)"
```

Run performance tests:

```bash
uv run pytest backend/fpbase/tests/test_performance.py -v
```

#### Step 4.11: Manual Testing Checklist (Enhanced)

**Frontend Dropdowns**:

- [ ] Microscope form detector dropdown opens/closes
- [ ] Microscope form allows clear selection
- [ ] FRET calculator donor selection works
- [ ] FRET calculator acceptor selection works
- [ ] Dropdowns match Bootstrap 4 theme styling
- [ ] Keyboard navigation works (Tab, Arrow keys, Enter, Esc)
- [ ] Clear button appears only when item selected
- [ ] Placeholder text displays when empty

**Django Forms**:

- [ ] Lineage parent autocomplete searches correctly
- [ ] Lineage parent autocomplete displays results
- [ ] Lineage parent autocomplete selects value
- [ ] Spectrum owner_state autocomplete works
- [ ] Optical config filter multi-select works
- [ ] Can add multiple filters
- [ ] Can remove individual filters (remove button on each item)
- [ ] Clear button works on single-select fields
- [ ] Autocomplete respects minimum_query_length (2 chars)
- [ ] Loading indicator appears during AJAX requests
- [ ] No results message displays appropriately
- [ ] Form pre-population works for edit forms
- [ ] Initial value displays correctly on page load

**Form Validation**:

- [ ] All forms submit successfully
- [ ] Selected values persist after submission
- [ ] Validation errors display correctly
- [ ] `.is-invalid` class applied to Tom-Select wrapper on error
- [ ] `.is-valid` class applied on success
- [ ] Error messages positioned correctly
- [ ] Form errors don't break Tom-Select widgets

**Accessibility** (WCAG 2.1 AA):

- [ ] ARIA labels present and correct (`aria-label`, `aria-labelledby`)
- [ ] ARIA roles assigned (`role="combobox"`, `role="listbox"`, `role="option"`)
- [ ] ARIA states update (`aria-expanded`, `aria-selected`, `aria-activedescendant`)
- [ ] Screen reader announces dropdown open/close
- [ ] Screen reader announces selected item changes
- [ ] Focus indicator visible (3:1 contrast minimum)
- [ ] Focus trap works when dropdown open
- [ ] Focus returns to input after selection
- [ ] High contrast mode displays correctly
- [ ] Touch targets ≥44x44px (mobile)
- [ ] Passes automated axe-core scan

**Performance**:

- [ ] Bundle sizes reduced (check webpack output)
- [ ] Page load times maintained or improved
- [ ] Autocomplete responsive (<200ms p95)
- [ ] No console errors
- [ ] No memory leaks (repeated interactions)
- [ ] Large datasets render smoothly (>100 options)

**Browser Compatibility**:

- [ ] Chrome (latest)
- [ ] Firefox (latest)
- [ ] Safari (latest)
- [ ] Edge (latest)
- [ ] Mobile Safari (iOS 15+)
- [ ] Chrome Android (latest)

---

## Deployment Strategy

### Pre-Deployment Checklist

**Code Ready**:

- [ ] All phases (0-4) complete
- [ ] All tests passing (unit, integration, E2E, a11y, performance)
- [ ] Code coverage ≥95%
- [ ] Bundle size reduced ≥15KB (verified)
- [ ] Manual testing checklist 100% complete
- [ ] Cross-browser testing passed
- [ ] Performance baselines met or exceeded
- [ ] Sentry alerts configured

**Staging Verification**:

- [ ] Deploy to staging environment
- [ ] Run full test suite in staging
- [ ] Manual smoke test in staging
- [ ] Test rollback procedure in staging (CRITICAL)
- [ ] Performance verification in staging
- [ ] Accessibility scan in staging

### Deployment Approach: Incremental Rollout with Feature Flags

**Strategy**: Gradual rollout with instant rollback capability via `USE_TOMSELECT` environment variable.

#### Phase A: Deploy Code (Feature Disabled)

```bash
# Deploy to production with feature flag OFF
heroku config:set USE_TOMSELECT=false

# Deploy code
git push heroku main

# Verify deployment successful
heroku ps
heroku logs --tail
```

At this point:
- Tom-Select code is deployed but not active
- Select2 still in use (zero user impact)
- Rollback = no-op (nothing changed)

#### Phase B: Enable for Staff Only (10% Traffic)

```python
# backend/config/settings/production.py

USE_TOMSELECT = env.bool("USE_TOMSELECT", default=False)

# OR percentage-based rollout
import random
USE_TOMSELECT = env.bool("USE_TOMSELECT", default=False) or (
    random.random() < 0.10  # 10% of requests
)

# OR staff-only
from django.contrib.auth import get_user
USE_TOMSELECT = lambda request: (
    env.bool("USE_TOMSELECT", default=False) or
    (hasattr(request, 'user') and request.user.is_staff)
)
```

Enable for staff:

```bash
heroku config:set USE_TOMSELECT=true

# Monitor Sentry for errors tagged with tomselect_migration:active
# Monitor for 24 hours
```

**Success Criteria** (24 hours):
- Zero critical errors in Sentry
- No user complaints
- Performance within ±5% baseline
- Bundle size reduced as expected

**Rollback Trigger**:
- Error rate >2x baseline → Immediate rollback
- Critical form submission failures → Immediate rollback
- >5 user complaints in 1 hour → Immediate rollback

#### Phase C: 50% Rollout

If Phase B successful:

```bash
# Increase to 50% of traffic
# (Update percentage-based logic or use feature flag service)

# Monitor for 48 hours
```

#### Phase D: 100% Rollout

If Phase C successful:

```bash
heroku config:set USE_TOMSELECT=true

# Monitor for 7 days before declaring success
```

#### Phase E: Remove Feature Flag (After 14 Days Stable)

If no issues for 14 days:

1. Remove feature flag from settings
2. Remove old Select2/DAL code (already done in Phase 4)
3. Deploy final cleanup

### Monitoring During Rollout

**Sentry Dashboard**:

Create custom Sentry dashboard:
- Total errors (grouped by tomselect_migration tag)
- Error rate trend
- Most common errors
- Affected users

**Metrics to Watch**:

1. **Error Rate**: Should stay ≤ baseline + 10%
2. **Form Submission Success Rate**: Should stay ≥99%
3. **Page Load Time**: Should improve 5-10%
4. **Bundle Size**: Should decrease ~19KB
5. **User Complaints**: Should be zero

**Alert Thresholds**:

| Metric | Warning | Critical | Action |
|--------|---------|----------|--------|
| Error Rate | +20% | +50% | Rollback immediately |
| Form Failures | >1% | >5% | Investigate / Rollback |
| Page Load Degradation | +10% | +25% | Investigate |
| Critical Bugs | 1 | 3 | Fix immediately or rollback |

### Rollback Procedure

**Instant Rollback** (if USE_TOMSELECT env var still exists):

```bash
# Disable feature flag (takes effect immediately)
heroku config:set USE_TOMSELECT=false

# Verify rollback
heroku logs --tail | grep "tomselect"
```

**Full Rollback** (if feature flag removed):

```bash
# Revert to previous release
heroku releases
heroku rollback v123  # Replace with last known good version

# Verify rollback
heroku ps
heroku logs --tail
```

**Post-Rollback**:

1. Analyze Sentry errors
2. Identify root cause
3. Fix in feature branch
4. Re-test in staging
5. Re-attempt deployment

### Post-Deployment

**Success Declaration** (after 14 days stable):

1. ✅ Zero critical errors in Sentry
2. ✅ No user complaints
3. ✅ Performance metrics met
4. ✅ Bundle size reduced as expected
5. ✅ WCAG 2.1 AA compliance maintained
6. ✅ All forms functioning correctly

**Document**:

- Update ADR status to "Implemented"
- Remove embedscope.js TODO from backlog (or schedule fix)
- Plan Bootstrap 5 migration (ADR-002)
- Share success metrics with team

---

## Code Examples

### Example 1: Basic Single Select

```javascript
import TomSelect from 'tom-select'

new TomSelect('#basic-select', {
  placeholder: 'Select an option...',
  plugins: ['clear_button'],
})
```

### Example 2: Multi-Select with Remove Buttons

```javascript
new TomSelect('#multi-select', {
  maxItems: null,  // Unlimited
  plugins: ['remove_button'],
  placeholder: 'Select multiple...',
})
```

### Example 3: AJAX Autocomplete with Error Handling

```javascript
new TomSelect('#ajax-select', {
  valueField: 'id',
  labelField: 'name',
  searchField: 'name',

  load: function(query, callback) {
    if (!query.length) return callback()

    fetch(`/api/search?q=${encodeURIComponent(query)}`)
      .then(response => {
        if (!response.ok) throw new Error(`HTTP ${response.status}`)
        return response.json()
      })
      .then(json => callback(json.results || []))
      .catch(error => {
        console.error('Autocomplete error:', error)
        callback()  // Empty results on error
      })
  },
})
```

### Example 4: Django Autocomplete View

```python
from django_tomselect.autocompletes import AutocompleteModelView

class ProteinAutocompleteView(AutocompleteModelView):
    model = Protein
    search_lookups = ["name__icontains", "slug__icontains"]
    ordering = ["name"]
    page_size = 20

    def create_result_dict(self, result):
        """Customize id/text fields."""
        return {
            "id": result.slug,
            "text": f"{result.name} ({result.slug})",
        }

    def get_queryset(self):
        """Apply custom filtering."""
        qs = super().get_queryset()
        if self.request.GET.get("fluorescent_only"):
            qs = qs.filter(is_fluorescent=True)
        return qs
```

### Example 5: Django Form Field

```python
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginClearButton

class MyForm(forms.ModelForm):
    protein = TomSelectModelChoiceField(
        config=TomSelectConfig(
            url="protein-autocomplete",
            placeholder="Select a protein...",
            minimum_query_length=2,
            highlight=True,
            preload="focus",
            plugin_clear_button=PluginClearButton(title="Clear"),
            css_framework="bootstrap4",
        ),
        queryset=Protein.objects.all(),
        required=False,
    )
```

---

## Testing Strategy

### Test-Driven Development Approach

**Principle**: Write FAILING tests BEFORE code changes (red-green-refactor).

**Example TDD Cycle for Phase 1**:

```python
# Step 1: Write failing test (RED)
@pytest.mark.xfail(reason="Tom-Select not implemented yet")
def test_microscope_form_uses_tomselect(page, live_server):
    page.goto(f'{live_server.url}/microscopes/add/')
    assert page.query_selector('.ts-wrapper') is not None

# Run: Should FAIL (xfail status)
uv run pytest -k test_microscope_form_uses_tomselect

# Step 2: Implement Tom-Select in microscope-form.js (GREEN)
# ... code changes ...

# Step 3: Remove xfail, run test: Should PASS
def test_microscope_form_uses_tomselect(page, live_server):  # No xfail
    page.goto(f'{live_server.url}/microscopes/add/')
    assert page.query_selector('.ts-wrapper') is not None

# Run: Should PASS
uv run pytest -k test_microscope_form_uses_tomselect

# Step 4: Refactor for quality, tests still pass
```

### Unit Tests

```python
# backend/proteins/tests/test_forms.py

class TestLineageFormMigration:
    def test_parent_field_uses_tomselect(self):
        """Verify parent field uses TomSelectModelChoiceField."""
        form = LineageForm()
        from django_tomselect.forms import TomSelectModelChoiceField
        assert isinstance(form.fields['parent'], TomSelectModelChoiceField)

    def test_form_submission_with_parent(self):
        """Verify form still submits correctly."""
        protein = Protein.objects.create(name="EGFP")
        parent_protein = Protein.objects.create(name="GFP")
        parent_lineage = Lineage.objects.create(protein=parent_protein)

        data = {
            'protein': protein.id,
            'parent': parent_lineage.id,
            'mutation': 'F64L'
        }

        form = LineageForm(data=data)
        assert form.is_valid()
        lineage = form.save()
        assert lineage.parent == parent_lineage
```

### Integration Tests

```python
# backend/proteins/tests/test_autocomplete.py

class TestProteinTomSelectView:
    def test_autocomplete_search(self, client):
        """Test search returns matching proteins."""
        Protein.objects.create(name="GFP", slug="gfp")

        response = client.get('/proteins/ts/protein/', {'q': 'GFP'})
        assert response.status_code == 200

        data = response.json()
        assert 'results' in data
        assert len(data['results']) == 1
        assert data['results'][0]['text'] == 'GFP'
        assert data['results'][0]['id'] == 'gfp'  # Slug, not pk
```

### End-to-End Tests

```python
# backend/fpbase/tests/test_end2end.py

def test_microscope_form_tom_select(page, live_server):
    """Test microscope form with Tom-Select."""
    page.goto(f'{live_server.url}/microscopes/add/')

    # Wait for Tom-Select
    detector_control = page.wait_for_selector('.ts-wrapper .ts-control')

    # Click to open
    detector_control.click()

    # Wait for dropdown
    dropdown = page.wait_for_selector('.ts-dropdown', state='visible')

    # Select first option
    first_option = dropdown.query_selector('.option')
    first_option.click()

    # Verify selection
    selected_item = detector_control.query_selector('.item')
    assert selected_item is not None
```

### Accessibility Tests

```python
# backend/fpbase/tests/test_accessibility.py

from axe_playwright_python import Axe

class TestAccessibility:
    def test_wcag_aa_compliance(self, page, live_server):
        """All forms must pass WCAG 2.1 AA."""
        page.goto(f'{live_server.url}/microscopes/add/')

        axe = Axe(page)
        results = axe.run()

        violations = results['violations']
        assert len(violations) == 0, f"WCAG violations: {axe.report(violations)}"

    def test_keyboard_navigation(self, page, live_server):
        """Complete form fillout using only keyboard."""
        page.goto(f'{live_server.url}/proteins/lineage/add/')

        # Tab to parent field
        page.keyboard.press('Tab')
        page.keyboard.press('Tab')

        # Open with Enter
        page.keyboard.press('Enter')

        # Navigate with arrows
        page.keyboard.press('ArrowDown')
        page.keyboard.press('ArrowDown')

        # Select with Enter
        page.keyboard.press('Enter')

        # Verify selection
        selected = page.query_selector('.ts-wrapper .item')
        assert selected is not None
```

### Performance Tests

```python
# backend/fpbase/tests/test_performance.py

import time

@pytest.mark.performance
class TestPerformance:
    def test_dropdown_open_latency(self, page, live_server):
        """Dropdown should open in <50ms."""
        page.goto(f'{live_server.url}/microscopes/add/')

        start = time.perf_counter()
        page.click('.ts-control')
        page.wait_for_selector('.ts-dropdown[style*="display"]')
        duration = (time.perf_counter() - start) * 1000

        assert duration < 50
```

---

## Rollback Plan

### Per-Phase Rollback

Each phase can be rolled back independently.

#### Phase 1 Rollback

```bash
# Revert microscope-form.js and fret.js
git checkout main -- frontend/src/microscope-form.js frontend/src/js/fret.js

# Reinstall Select2 if removed
pnpm add select2 select2-theme-bootstrap4

# Rebuild
pnpm build
```

#### Phase 2 Rollback (Feature Flag)

```bash
# Instant rollback via environment variable
heroku config:set USE_TOMSELECT=false

# OR full code rollback
git checkout main -- backend/proteins/forms/ backend/proteins/views/
uv remove django-tomselect
uv add django-autocomplete-light
uv sync
```

#### Phase 3 Rollback

```bash
git checkout main -- frontend/src/js/project.js
pnpm build
```

#### Phase 4 Rollback (Full Revert)

```bash
# Revert all changes
git checkout main -- frontend/ backend/

# Restore dependencies
pnpm install
uv sync

# Rebuild
pnpm build
uv run backend/manage.py collectstatic --noinput
```

### Rollback Testing

**CRITICAL**: Test rollback in staging BEFORE production deployment:

```bash
# In staging environment:
# 1. Deploy migration
# 2. Verify working
# 3. Execute rollback procedure
# 4. Verify rollback successful
# 5. Re-deploy migration
# 6. Document any issues
```

---

## Appendices

### Appendix A: Select2 vs Tom-Select API Comparison

| Feature | Select2 | Tom-Select |
|---------|---------|------------|
| Initialization | `$('#select').select2({...})` | `new TomSelect('#select', {...})` |
| AJAX | `ajax: { url: '...' }` | `load: function(query, callback) {...}` |
| Placeholder | `placeholder: 'text'` | `placeholder: 'text'` |
| Clear button | `allowClear: true` | `plugins: ['clear_button']` |
| Multi-select | `multiple: true` (HTML) | `maxItems: null` |
| Theming | `theme: 'bootstrap'` | Import bootstrap CSS |
| Value field | `data-value-field` | `valueField: 'id'` |
| Label field | `data-text-field` | `labelField: 'name'` |
| Dependencies | Requires jQuery | Zero dependencies |

### Appendix B: django-autocomplete-light vs django-tomselect

| Feature | DAL | django-tomselect |
|---------|-----|------------------|
| View base | `autocomplete.Select2QuerySetView` | `AutocompleteModelView` |
| Form field | `autocomplete.ModelSelect2` widget | `TomSelectModelChoiceField` |
| Multi-select | `ModelSelect2Multiple` | `TomSelectModelMultipleChoiceField` |
| Configuration | Widget attrs | `TomSelectConfig` object |
| Result customization | `get_results()` | `create_result_dict()` ✅ |
| URL in widget | `url='namespace:name'` | `url='namespace:name'` (in config) |
| Theming | `data-theme='bootstrap'` | `css_framework='bootstrap4'` |
| Dependencies | Select2 + jQuery | Tom-Select (no jQuery) |

### Appendix C: Bundle Size Analysis (CORRECTED)

**Before Migration**:

| Component | Size (gzipped) |
|-----------|----------------|
| Select2 | 30KB |
| jQuery | 30KB |
| select2-theme-bootstrap4 | 5KB |
| Custom _select2-bootstrap.scss | 5KB |
| **Total** | **70KB** |

**After Migration**:

| Component | Size (gzipped) |
|-----------|----------------|
| Tom-Select | 16KB |
| jQuery | 30KB ⚠️ (STAYS - Bootstrap 4 requires it) |
| Tom-Select Bootstrap CSS | 3KB |
| Custom validation styles | 2KB |
| **Total** | **51KB** |

**Net Savings**: ~19KB per bundle (~27% reduction, NOT 73%)

**Important**: jQuery stays at 30KB because:
- Bootstrap 4 JavaScript (modals, tooltips, etc.) requires jQuery
- `jquery.formset.js` requires jQuery
- 375+ jQuery calls in frontend/src/js/

**Full jQuery removal** requires Bootstrap 5 migration (future ADR-002).

### Appendix D: Sentry Issue FPBASE-5H3 Root Cause

**Error**: `TypeError: Cannot read properties of undefined (reading 'toString')`

**Location**: `select2@4.0.13/dist/js/select2.full.js:965:21`

**Code**:

```javascript
var selectedIds = $.map(selected, function (s) {
  return s.id.toString();  // ❌ Crashes when s.id is undefined
});
```

**Root Cause**: Select2 expects all selected items to have `id` property, but encounters item where `s.id` is `undefined`.

**Why**: Race condition during initialization or malformed server response.

**How Tom-Select Fixes This**:

- Better error handling
- More defensive code for edge cases
- Modern, actively maintained codebase
- No reliance on jQuery's `$.map` behavior

---

## Resources

### Documentation

**Tom-Select**:
- Website: https://tom-select.js.org/
- GitHub: https://github.com/orchidjs/tom-select
- Docs: https://tom-select.js.org/docs/
- Examples: https://tom-select.js.org/examples/
- Plugins: https://tom-select.js.org/docs/plugins/

**django-tomselect**:
- PyPI: https://pypi.org/project/django-tomselect/
- GitHub: https://github.com/OmenApps/django-tomselect
- Docs: https://django-tomselect.readthedocs.io/
- Tutorial: https://django-tomselect.readthedocs.io/en/latest/tutorial.html
- Config: https://django-tomselect.readthedocs.io/en/latest/api/config.html

### Community

- Tom-Select Issues: https://github.com/orchidjs/tom-select/issues
- django-tomselect Issues: https://github.com/OmenApps/django-tomselect/issues

---

## Migration Completion Checklist

**Before Starting**:

- [ ] Team briefed on migration plan
- [ ] All expert reviews incorporated
- [ ] Feature branch created: `feature/select2-to-tomselect`
- [ ] Staging environment ready

**Phase 0 Complete**:

- [ ] Performance baselines collected
- [ ] Accessibility testing setup (axe-core installed)
- [ ] TDD test stubs created (failing tests)
- [ ] Feature flag implemented
- [ ] Sentry alerts configured
- [ ] Cross-browser testing configured

**Phase 1 Complete**:

- [ ] Tom-Select installed
- [ ] microscope-form.js migrated and tested
- [ ] fret.js migrated and tested
- [ ] embedscope.js strategy documented (keeps Select2)
- [ ] Bundle size reduced (verified)
- [ ] TDD tests passing

**Phase 2 Complete**:

- [ ] django-tomselect installed
- [ ] Settings configured (no middleware/context processor!)
- [ ] Autocomplete views created with `create_result_dict()`
- [ ] URL patterns added (parallel /ts/* endpoints)
- [ ] LineageForm migrated and tested
- [ ] SpectrumForm migrated and tested
- [ ] OpticalConfigForm migrated and tested
- [ ] ReferenceAutocomplete migrated and tested
- [ ] All form tests passing

**Phase 3 Complete**:

- [ ] project.js migrated with error handling
- [ ] AJAX autocomplete tested

**Phase 4 Complete**:

- [ ] Select2/theme removed from package.json
- [ ] django-autocomplete-light removed
- [ ] Old DAL views/URLs deleted
- [ ] SCSS cleaned up (validation styles kept)
- [ ] E2E tests updated (all Select2 selectors migrated)
- [ ] Full test suite passing (≥95% coverage)
- [ ] Bundle sizes verified (<baseline - 15KB)
- [ ] Accessibility scan clean (zero WCAG violations)
- [ ] Cross-browser tests passing
- [ ] Performance tests passing

**Deployment Complete**:

- [ ] Staging deployment successful
- [ ] Rollback tested in staging
- [ ] Production deployment (feature disabled)
- [ ] Staff rollout successful (24 hours)
- [ ] 50% rollout successful (48 hours)
- [ ] 100% rollout successful (7 days)
- [ ] Feature flag removed (14 days stable)
- [ ] ADR updated to "Implemented"
- [ ] Success metrics documented
- [ ] Team celebration! 🎉

---

**End of Migration Plan**

This corrected guide incorporates all expert feedback and is ready for immediate implementation by an AI code agent or development team.
