# Select2 → Tom-Select Migration: Implementation Guide

**Status**: Ready for Execution
**Last Updated**: 2025-01-26

## Phase 0: Pre-Migration Setup

### Step 0.1: Create Baseline Metrics

```bash
# Create baseline file
echo "=== Baseline Metrics $(date) ===" > baseline_metrics.txt

# Bundle sizes
pnpm build
echo "Bundle sizes:" >> baseline_metrics.txt
ls -lh backend/fpbase/static/bundles/*.js >> baseline_metrics.txt

# Current test status
uv run pytest --co -q > test_inventory.txt
echo "Test count: $(wc -l < test_inventory.txt)" >> baseline_metrics.txt
```

### Step 0.2: Install Accessibility Testing

```bash
pnpm add --save-dev @axe-core/playwright axe-core
```

Create test file:

```python
# backend/fpbase/tests/test_accessibility.py

import pytest
from axe_playwright_python import Axe

@pytest.fixture
def axe(page):
    return Axe(page)

@pytest.mark.accessibility
class TestAccessibility:
    def test_microscope_form_wcag(self, page, live_server, axe):
        page.goto(f'{live_server.url}/microscopes/add/')
        results = axe.run()
        assert len(results['violations']) == 0, f"WCAG violations: {axe.report(results['violations'])}"
```

Run baseline:

```bash
uv run pytest backend/fpbase/tests/test_accessibility.py -v
# Document any existing violations
```

### Step 0.3: Create TDD Test Stubs

```python
# backend/fpbase/tests/test_tomselect_migration.py

import pytest

class TestTomSelectMigration:
    @pytest.mark.xfail(reason="Tom-Select not yet implemented")
    def test_microscope_form_uses_tomselect(self, page, live_server):
        page.goto(f'{live_server.url}/microscopes/add/')
        assert page.query_selector('.ts-wrapper') is not None
        assert page.query_selector('.select2-container') is None

    @pytest.mark.xfail(reason="django-tomselect not yet implemented")
    def test_lineage_form_uses_tomselect_field(self):
        from proteins.forms.forms import LineageForm
        from django_tomselect.forms import TomSelectModelChoiceField
        form = LineageForm()
        assert isinstance(form.fields['parent'], TomSelectModelChoiceField)
```

Verify they fail:

```bash
uv run pytest backend/fpbase/tests/test_tomselect_migration.py -v
# Should show 2 xfail
```

### Step 0.4: Add Feature Flag

```python
# backend/config/settings/base.py

# Add after other env variables
USE_TOMSELECT = env.bool("USE_TOMSELECT", default=False)
```

### Step 0.5: Configure Sentry Alert (Production Only)

If using Sentry, add to production config:

```python
# backend/config/settings/production.py

# In sentry_sdk.init before_send callback, add:
"tags": {
    **event.get("tags", {}),
    "tomselect_migration": "active" if USE_TOMSELECT else "inactive",
}
```

Set up alert in Sentry UI:

- Filter: `tags.tomselect_migration:active AND message contains "TomSelect"`
- Threshold: >10 events/hour
- Action: Alert team

---

## Phase 1: Simple Frontend Dropdowns

### Step 1.1: Install Tom-Select

```bash
pnpm add tom-select
```

### Step 1.2: Update microscope-form.js

```javascript
// frontend/src/microscope-form.js

// REMOVE these imports:
// import "select2/dist/css/select2.css"
// import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
// import "select2/dist/js/select2.full"

// ADD these imports:
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

// REPLACE this code:
// $("#id_detector, #id_light_source, #id_extra_cameras, #id_extra_lights").select2({
//   theme: "bootstrap",
//   containerCssClass: ":all:",
//   placeholder: "---------",
//   allowClear: true,
//   width: "auto",
// })

// WITH this code:
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
    controlClass: element.className,
  })
})
```

### Step 1.3: Update fret.js

```javascript
// frontend/src/js/fret.js

// REMOVE Select2 initialization (around lines 113-122)

// ADD at appropriate location:
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'

['#donor-select', '#acceptor-select'].forEach(selector => {
  const element = document.querySelector(selector)
  if (!element) return

  new TomSelect(selector, {
    controlClass: element.className,
  })
})
```

### Step 1.4: Update embedscope.js (KEEP Select2)

```javascript
// frontend/src/embedscope.js

// KEEP these imports (DO NOT REMOVE - legacy microscope.js needs Select2):
import "select2/dist/js/select2.full.js"

// ADD these imports for new features:
import TomSelect from 'tom-select'
import 'tom-select/dist/css/tom-select.bootstrap4.css'
```

Document tech debt:

```markdown
<!-- Add to project docs/technical-debt.md -->
## embedscope.js Select2 Dependency

**Issue**: embedscope.js must keep Select2 because legacy microscope.js (CDN-loaded) depends on it
**Blocker**: backend/fpbase/static/js/microscope.js (847 lines)
**To Fix**: Refactor microscope.js to use Tom-Select OR bundle with webpack
**Effort**: 2-3 days
**Priority**: Low
```

### Step 1.5: Build and Test

```bash
# Build
pnpm build

# Run TDD test (remove xfail first)
# Edit backend/fpbase/tests/test_tomselect_migration.py:
# Remove @pytest.mark.xfail decorator from test_microscope_form_uses_tomselect

uv run pytest backend/fpbase/tests/test_tomselect_migration.py::TestTomSelectMigration::test_microscope_form_uses_tomselect -v
# Should PASS

# Run E2E tests
uv run pytest backend/fpbase/tests/test_end2end.py::TestPagesRender::test_microscopes -v
uv run pytest backend/fpbase/tests/test_end2end.py::TestPagesRender::test_fret -v
```

### Step 1.6: Manual Verification

Start dev server:

```bash
pnpm dev
```

Test:

1. Go to `/microscopes/add/`
2. Click detector dropdown → opens
3. Select option → works
4. Click X (clear) → clears
5. Go to `/fret/`
6. Test donor/acceptor dropdowns

### Step 1.7: Check Bundle Size

```bash
pnpm build
ls -lh backend/fpbase/static/bundles/microscope-form.*.js
ls -lh backend/fpbase/static/bundles/fret.*.js
# Compare to baseline_metrics.txt
```

---

## Phase 2: Django Integration

### Step 2.1: Install django-tomselect

```bash
uv add django-tomselect
```

### Step 2.2: Update Settings

```python
# backend/config/settings/base.py

INSTALLED_APPS = [
    # ... existing apps
    "dal",  # Keep during migration
    "dal_select2",
    "django_tomselect",  # ADD this
]

# Optional global config:
TOMSELECT_CONFIG = {
    "css_framework": "bootstrap4",
}

# DO NOT ADD middleware or context processors - they don't exist!
```

### Step 2.3: Collect Static

```bash
uv run backend/manage.py collectstatic --noinput
```

### Step 2.4: Create Autocomplete Views

```python
# backend/proteins/views/autocomplete.py

# ADD these imports at top:
from __future__ import annotations
from django.contrib.auth.mixins import LoginRequiredMixin
from django_tomselect.autocompletes import AutocompleteModelView

# ADD these NEW view classes (keep old DAL views for now):

class ProteinTomSelectView(AutocompleteModelView):
    model = Protein
    search_lookups = ["name__icontains", "slug__icontains"]
    ordering = ["name"]
    page_size = 20

    def create_result_dict(self, result):
        return {
            "id": result.slug,
            "text": result.name,
        }

    def get_queryset(self):
        qs = super().get_queryset()
        if self.request.GET.get("type") == "spectra":
            qs = qs.filter(default_state__spectra__isnull=False).distinct()
        return qs


class LineageTomSelectView(AutocompleteModelView):
    model = Lineage
    search_lookups = ["protein__name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        return {
            "id": result.pk,
            "text": result.protein.name,
        }


class StateTomSelectView(AutocompleteModelView):
    model = State
    search_lookups = ["protein__name__icontains", "name__icontains"]
    ordering = ["protein__name"]
    page_size = 20

    def get_queryset(self):
        return super().get_queryset().select_related("protein")

    def create_result_dict(self, result):
        return {
            "id": result.pk,
            "text": str(result),
        }


class FilterTomSelectView(AutocompleteModelView):
    model = Filter
    search_lookups = ["name__icontains", "part__icontains"]
    ordering = ["name"]
    page_size = 20
```

### Step 2.5: Add URL Patterns

```python
# backend/proteins/urls.py

# ADD to imports:
from .views.autocomplete import (
    # ... existing DAL imports
    ProteinTomSelectView,
    LineageTomSelectView,
    StateTomSelectView,
    FilterTomSelectView,
)

# ADD to urlpatterns (keep old DAL URLs for now):
urlpatterns = [
    # ... existing patterns

    # New Tom-Select endpoints
    path("ts/protein/", ProteinTomSelectView.as_view(), name="protein-ts"),
    path("ts/lineage/", LineageTomSelectView.as_view(), name="lineage-ts"),
    path("ts/state/", StateTomSelectView.as_view(), name="state-ts"),
    path("ts/filter/", FilterTomSelectView.as_view(), name="filter-ts"),
]
```

### Step 2.6: Migrate LineageForm

```python
# backend/proteins/forms/forms.py

# ADD to imports:
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginClearButton

# FIND the LineageForm class and REPLACE parent field:

class LineageForm(forms.ModelForm):
    # REPLACE this:
    # parent = forms.ModelChoiceField(
    #     queryset=Lineage.objects.all().prefetch_related("protein"),
    #     widget=autocomplete.ModelSelect2(...),
    #     ...
    # )

    # WITH this:
    parent = TomSelectModelChoiceField(
        config=TomSelectConfig(
            url="proteins:lineage-ts",
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

    # Keep existing Meta, clean(), __init__() methods unchanged
```

Test:

```bash
# Remove xfail from test
# Edit backend/fpbase/tests/test_tomselect_migration.py
# Remove @pytest.mark.xfail from test_lineage_form_uses_tomselect_field

uv run pytest backend/fpbase/tests/test_tomselect_migration.py::TestTomSelectMigration::test_lineage_form_uses_tomselect_field -v
# Should PASS

uv run pytest backend/proteins/tests/test_forms.py::TestLineageForm -v
# Should PASS
```

### Step 2.7: Migrate SpectrumForm

```python
# backend/proteins/forms/spectrum.py

# ADD to imports:
from django_tomselect.forms import TomSelectModelChoiceField
from django_tomselect.app_settings import TomSelectConfig

# FIND SpectrumForm and REPLACE owner_state field:

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

Test:

```bash
uv run pytest backend/proteins/tests/test_forms.py -k Spectrum -v
```

### Step 2.8: Migrate OpticalConfigForm

```python
# backend/proteins/forms/microscope.py

# ADD to imports:
from django_tomselect.forms import TomSelectModelMultipleChoiceField
from django_tomselect.app_settings import TomSelectConfig, PluginRemoveButton

# FIND MultipleFilterField and REPLACE:

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
                max_items=None,
                plugin_remove_button=PluginRemoveButton(title="Remove"),
                css_framework="bootstrap4",
            ),
            queryset=Filter.objects.all(),
            required=False,
        )

# OpticalConfigForm should use MultipleFilterField as before (no changes needed)
```

Test:

```bash
uv run pytest backend/proteins/tests/test_forms.py -k OpticalConfig -v
```

### Step 2.9: Migrate ReferenceAutocomplete

```python
# backend/references/views.py

# ADD to imports:
from django.contrib.auth.mixins import LoginRequiredMixin
from django_tomselect.autocompletes import AutocompleteModelView

# ADD new view:

class ReferenceTomSelectView(LoginRequiredMixin, AutocompleteModelView):
    model = Reference
    search_lookups = ["doi__icontains"]
    ordering = ["citation"]
    page_size = 20

    def create_result_dict(self, result):
        return {
            "id": result.doi,
            "text": result.citation,
        }
```

Add URL in references/urls.py:

```python
path("ts/reference/", ReferenceTomSelectView.as_view(), name="reference-ts"),
```

Update reference form to use Tom-Select field.

### Step 2.10: Verify Template Includes form.media

Check templates:

```bash
grep -r "{{ form.media }}" backend/fpbase/templates/
```

If missing, add to form templates:

```django
{% block extra_css %}
    {{ form.media.css }}
{% endblock %}

{% block extra_js %}
    {{ form.media.js }}
{% endblock %}
```

### Step 2.11: Run Full Form Tests

```bash
uv run pytest backend/proteins/tests/test_forms.py -v
uv run pytest backend/references/tests/ -v
```

---

## Phase 3: Custom AJAX

### Step 3.1: Update project.js

```javascript
// frontend/src/js/project.js

// FIND Select2 initialization (around line 207-228)
// REPLACE with:

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
          throw new Error(`HTTP ${response.status}`)
        }
        return response.json()
      })
      .then(json => {
        callback(json.results || [])
      })
      .catch(error => {
        console.error('Autocomplete error:', error)
        callback()
      })
  },

  render: {
    option: function(item, escape) {
      return `<div>${escape(item.text)}</div>`
    }
  }
})
```

### Step 3.2: Test

```bash
pnpm build

# Manual test:
# Navigate to protein slug builder page
# Type in autocomplete
# Verify results
# Select protein
# Verify URL building works
```

---

## Phase 4: Cleanup & Final Testing

### Step 4.1: Remove Old Dependencies

```bash
# Frontend (DO NOT remove jQuery!)
pnpm remove select2 select2-theme-bootstrap4

# Backend
uv remove django-autocomplete-light
uv sync
```

### Step 4.2: Delete Old Code

Delete old DAL views:

```python
# In backend/proteins/views/autocomplete.py
# DELETE these classes:
# - ProteinAutocomplete
# - LineageAutocomplete
# - StateAutocomplete
# - FilterAutocomplete

# DELETE old DAL imports if no longer used
```

Delete old URL patterns:

```python
# In backend/proteins/urls.py
# DELETE these paths:
# path("autocomplete-protein/", ...)
# path("lineage-autocomplete/", ...)
# path("state-autocomplete/", ...)
# path("filter-autocomplete/", ...)
```

Remove Select2 imports from JavaScript:

```bash
# Find Select2 imports (except embedscope.js)
grep -r "import.*select2" frontend/src/ --include="*.js" | grep -v embedscope

# Remove them manually (except in embedscope.js)
```

### Step 4.3: Update SCSS

```scss
// frontend/src/css/_tomselect-django.scss

// KEEP only validation state styles:
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

// DELETE all other Select2-Bootstrap styles
```

Rename and update imports:

```bash
# Rename file
mv frontend/src/css/_select2-bootstrap.scss frontend/src/css/_tomselect-django.scss

# Update main SCSS file imports
```

### Step 4.4: Update E2E Tests

```python
# backend/fpbase/tests/test_end2end.py

# FIND all Select2 selectors and REPLACE:

# .select2-container → .ts-wrapper
# .select2-selection → .ts-control
# .select2-search__field → .ts-control input[type="text"]
# #select2-{id}-container → #{id} + .ts-wrapper
# .select2-results → .ts-dropdown
# .select2-results__option → .ts-dropdown .option
# .select2-selection__clear → .clear-button
```

Example replacement:

```python
# Before:
# select2_container = page.query_selector('.select2-container')

# After:
ts_wrapper = page.query_selector('.ts-wrapper')

# Before:
# page.click('.select2-selection')

# After:
page.click('.ts-control')
```

### Step 4.5: Run Full Test Suite

```bash
uv run pytest --cov --cov-report=html --cov-fail-under=95
# Must PASS with ≥95% coverage
```

### Step 4.6: Accessibility Scan

```bash
uv run pytest backend/fpbase/tests/test_accessibility.py -v
# Must PASS with zero violations
```

### Step 4.7: Cross-Browser Testing

```bash
uv run pytest backend/fpbase/tests/test_end2end.py --browser chromium -v
uv run pytest backend/fpbase/tests/test_end2end.py --browser firefox -v
uv run pytest backend/fpbase/tests/test_end2end.py --browser webkit -v
# All must PASS
```

### Step 4.8: Build and Compare Bundle Sizes

```bash
pnpm build

echo "=== After Migration ===" >> bundle_comparison.txt
ls -lh backend/fpbase/static/bundles/*.js >> bundle_comparison.txt

# Compare with baseline_metrics.txt
diff baseline_metrics.txt bundle_comparison.txt
```

---

## Deployment

### Pre-Deployment Checklist

Run through checklist:

```bash
# All tests passing
uv run pytest

# Build successful
pnpm build

# No console errors in manual testing
# (Check browser console during manual tests)

# Migrations applied (if any)
uv run backend/manage.py migrate

# Static files collected
uv run backend/manage.py collectstatic --noinput
```

### Deploy to Staging

```bash
# Deploy code to staging (your deployment method)

# Verify in staging:
# - All forms work
# - Autocomplete works
# - Validation works
# - Styling correct
```

### Production Rollout

```bash
# Deploy with feature flag OFF
heroku config:set USE_TOMSELECT=false
git push heroku main

# Enable for staff only (24 hours)
heroku config:set USE_TOMSELECT=true

# Monitor Sentry for errors
# If issues: heroku config:set USE_TOMSELECT=false

# If stable, full rollout
# (Keep monitoring for 7 days)

# After 14 days stable, remove feature flag code
```

---

## Rollback Procedures

### Phase 1 Rollback

```bash
git checkout main -- frontend/src/microscope-form.js frontend/src/js/fret.js
pnpm add select2 select2-theme-bootstrap4
pnpm build
```

### Phase 2 Rollback (Feature Flag)

```bash
heroku config:set USE_TOMSELECT=false
# Instant rollback
```

### Phase 2 Rollback (Full)

```bash
git checkout main -- backend/proteins/forms/ backend/proteins/views/
uv remove django-tomselect
uv add django-autocomplete-light
uv sync
```

### Phase 3 Rollback

```bash
git checkout main -- frontend/src/js/project.js
pnpm build
```

### Full Rollback

```bash
git checkout main -- frontend/ backend/
pnpm install
uv sync
pnpm build
uv run backend/manage.py collectstatic --noinput
```

---

## Verification Commands

Quick verification commands to run at any point:

```bash
# Tests passing
uv run pytest -x

# Build works
pnpm build

# Server starts
uv run backend/manage.py check
uv run backend/manage.py runserver

# No import errors
uv run backend/manage.py shell -c "from proteins.forms.forms import LineageForm; print('OK')"

# Static files exist
ls backend/fpbase/static/bundles/*.js

# Coverage maintained
uv run pytest --cov --cov-report=term-missing
```

---

## Completion Checklist

**Phase 0**:

- [ ] Baseline metrics collected
- [ ] Accessibility tests added
- [ ] TDD stubs created (failing)
- [ ] Feature flag added

**Phase 1**:

- [ ] Tom-Select installed
- [ ] microscope-form.js updated
- [ ] fret.js updated
- [ ] embedscope.js documented (keeps Select2)
- [ ] TDD test passing
- [ ] Manual testing done

**Phase 2**:

- [ ] django-tomselect installed
- [ ] Settings updated (INSTALLED_APPS only)
- [ ] Autocomplete views created
- [ ] URLs added
- [ ] LineageForm migrated and tested
- [ ] SpectrumForm migrated and tested
- [ ] OpticalConfigForm migrated and tested
- [ ] ReferenceAutocomplete migrated and tested

**Phase 3**:

- [ ] project.js updated
- [ ] AJAX autocomplete tested

**Phase 4**:

- [ ] Old dependencies removed
- [ ] Old code deleted
- [ ] SCSS cleaned up
- [ ] E2E tests updated
- [ ] Full test suite passing (≥95% coverage)
- [ ] Accessibility scan clean
- [ ] Cross-browser tests passing

**Deployment**:

- [ ] Staging verified
- [ ] Production deployed (feature disabled)
- [ ] Staff rollout (24h stable)
- [ ] Full rollout (7d stable)
- [ ] Feature flag removed (14d stable)

---

**End of Implementation Guide**
