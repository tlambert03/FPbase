# FPbase Caching Overhaul Plan

**Date:** 2025-01-09
**Related Issue:** [#361 - Add ETag support for spectra list caching](https://github.com/tlambert03/FPbase/issues/361)

---

## Executive Summary

FPbase's current caching strategy is ad-hoc with arbitrary TTLs (60s, 1h, 6h, 24h), no systematic invalidation, and missing HTTP semantics (ETags, proper Vary headers). This document catalogs all existing caching patterns and proposes a comprehensive overhaul to a **signal-based, version-tracked system** with proper HTTP caching.

### Key Problems
- ❌ Arbitrary TTLs everywhere - no clear strategy
- ❌ No invalidation for most caches - data becomes stale
- ❌ Caching ORM objects directly (problematic)
- ❌ Manual cache invalidation requires remembering to call helpers
- ❌ Missing ETags - browsers re-download unchanged 140 KB payloads
- ❌ No Vary headers on API endpoints
- ❌ Unbounded cache growth in some areas

### Expected Outcomes
- ✅ **99.7% bandwidth reduction** for stable resources (140 KB → <1 KB validation)
- ✅ **Zero stale data** - automatic invalidation on DB changes
- ✅ **Zero manual invalidation** - signal-based automation
- ✅ **Better performance** - longer TTLs + ETags = fewer cache misses
- ✅ **Standards compliance** - proper HTTP semantics

---

## Table of Contents

1. [Current Caching Catalog](#current-caching-catalog)
2. [Critical Issues Summary](#critical-issues-summary)
3. [Implementation Plan](#implementation-plan)
4. [Design Decisions](#design-decisions)

---

## Current Caching Catalog

### 1. Django Cache Framework Usage

#### 1.1 Cache Backend Configuration

**Production (Redis):**
```python
# backend/config/settings/production.py:174-187
CACHES = {
    "default": {
        "BACKEND": "django_redis.cache.RedisCache",
        "LOCATION": REDIS_LOCATION,
        "OPTIONS": {
            "CLIENT_CLASS": "django_redis.client.DefaultClient",
            "IGNORE_EXCEPTIONS": True,  # mimics memcache behavior
            "CONNECTION_POOL_KWARGS": {"ssl_cert_reqs": ssl.CERT_NONE},
        },
    }
}
```

**Local (In-Memory):**
```python
# backend/config/settings/local.py:51-56
CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.locmem.LocMemCache",
        "LOCATION": "",
    }
}
```

---

#### 1.2 Direct Cache Operations (cache.get/set/delete)

##### A. GenBank Sequence Caching
**File:** `backend/proteins/extrest/entrez.py:124-145`

```python
def get_cached_gbseqs(gbids, max_age=60 * 60 * 24) -> dict[str, tuple[str, float]]:
    gbseqs = cache.get("gbseqs", {})
    now = time.time()
    tofetch = [id for id in gbids if id not in gbseqs or gbseqs[id][1] - now >= max_age]
    gbseqs.update({k: (v, now) for k, v in _fetch_gb_seqs(tofetch).items()})
    cache.set("gbseqs", gbseqs, 60 * 60 * 24)
    return gbseqs
```

- **What's cached:** GenBank protein sequences by ID
- **TTL:** 24 hours
- **Cache key:** `"gbseqs"`
- **Invalidation:** None (time-based expiry only)
- **Issues:**
  - ⚠️ Manual timestamp tracking
  - ⚠️ Cache grows unbounded (never removes old entries)
  - ⚠️ Arbitrary 24-hour TTL

---

##### B. Spectrum Data Caching (GraphQL)
**File:** `backend/proteins/schema/query.py:13-33`

```python
def get_cached_spectrum(id, timeout=60 * 60 * 24):
    key = f"_spectrum_{id}"
    spectrum = cache.get(key)
    if not spectrum:
        try:
            spectrum = (
                models.Spectrum.objects.filter(id=id)
                .select_related(...)
                .get()
            )
            cache.set(key, spectrum, timeout)
        except models.Spectrum.DoesNotExist:
            return None
    return spectrum
```

- **What's cached:** Individual spectrum ORM objects with related data
- **TTL:** 24 hours
- **Cache key:** `f"_spectrum_{id}"`
- **Invalidation:** None
- **Issues:**
  - 🔴 **CRITICAL:** Caching entire ORM objects can cause issues
  - ⚠️ No invalidation on spectrum updates

**Spectra List Caching:**
**File:** `backend/proteins/schema/query.py:98-122`

```python
def resolve_spectra(self, info, **kwargs):
    # ...
    if "owner" in requested_fields:
        key_suffix = "_".join(f"{k}_{v}" for k, v in sorted(fkwargs.items()))
        cache_key = f"_spectra_sluglist_{key_suffix}" if key_suffix else "_spectra_sluglist"

        result = cache.get(cache_key)
        if not result:
            result = models.Spectrum.objects.sluglist(filters=fkwargs or None)
            cache.set(cache_key, result, 60 * 60)  # 1 hour
        return result
```

- **What's cached:** Filtered spectrum sluglists
- **TTL:** 1 hour
- **Cache key:** `"_spectra_sluglist_{filters}"`
- **Invalidation:** None

---

##### C. Spectra Info Caching (REST API)
**File:** `backend/proteins/models/spectrum.py:72-80`

```python
SPECTRA_CACHE_KEY = "spectra_sluglist"

def get_cached_spectra_info(timeout=60 * 60):
    spectrainfo = cache.get(SPECTRA_CACHE_KEY)
    if not spectrainfo:
        spectrainfo = json.dumps({"data": {"spectra": Spectrum.objects.sluglist()}})
        cache.set(SPECTRA_CACHE_KEY, spectrainfo, timeout)
    return spectrainfo
```

- **What's cached:** JSON-serialized spectrum sluglist
- **TTL:** 1 hour
- **Cache key:** `"spectra_sluglist"`
- **Invalidation:** ✅ Called in `Spectrum.save()` (line 424) via `cache.delete(SPECTRA_CACHE_KEY)`
- **Used in:** `backend/proteins/api/views.py:31`
- **Issues:**
  - ⚠️ Manual cache deletion on save (error-prone)
  - ⚠️ Arbitrary 1-hour TTL
  - ⚠️ **Related to issue #361** - needs ETag support

---

##### D. Optical Config Caching
**File:** `backend/proteins/models/microscope.py:182-199`

```python
OC_CACHE_KEY = "optical_config_list"

def get_cached_optical_configs(timeout=60 * 60):
    ocinfo = cache.get(OC_CACHE_KEY)
    if not ocinfo:
        vals = OpticalConfig.objects.all().values(...)
        # ... processing ...
        ocinfo = json.dumps({"data": {"opticalConfigs": ocinfo}})
        cache.set(OC_CACHE_KEY, ocinfo, timeout)
    return ocinfo
```

- **What's cached:** JSON-serialized optical config data
- **TTL:** 1 hour
- **Cache key:** `"optical_config_list"`
- **Invalidation:** ✅ Called in `OpticalConfig.save()` (line 233) via `cache.delete(OC_CACHE_KEY)`
- **Used in:** `backend/proteins/api/views.py:36`
- **Issues:**
  - ⚠️ Manual cache deletion on save

---

##### E. Protein Slug Dictionary Caching
**File:** `backend/proteins/util/helpers.py:27-42`

```python
def link_excerpts(excerpts_qs, obj_name=None, aliases=()):
    # ...
    slug_dict = cache.get_or_set("slug_dict", create_slug_dict, 60)
```

- **What's cached:** Dictionary mapping protein names/aliases to slugs
- **TTL:** **60 seconds**
- **Cache key:** `"slug_dict"`
- **Invalidation:** None
- **Issues:**
  - 🔴 **CRITICAL:** 60-second TTL is way too short, causes frequent cache misses
  - ⚠️ No invalidation on protein updates

---

##### F. Google Analytics Popular Proteins
**File:** `backend/proteins/extrest/ga.py:44-55`

```python
def cached_ga_popular(max_age=60 * 60 * 24):
    results = cache.get("ga_popular_proteins")
    if not results:
        client = get_client()
        results = {
            "day": ga_popular_proteins(client, 1),
            "week": ga_popular_proteins(client, 7),
            "month": ga_popular_proteins(client, 30),
            "year": ga_popular_proteins(client, 365),
        }
        cache.set("ga_popular_proteins", results, max_age)
    return results
```

- **What's cached:** Google Analytics popular proteins data
- **TTL:** 24 hours
- **Cache key:** `"ga_popular_proteins"`
- **Invalidation:** None
- **Used in:** `backend/proteins/views/protein.py:524`
- **Issues:** None (appropriate for GA data - external data source)

---

##### G. FRET Calculation Caching
**File:** `backend/proteins/views/fret.py:17-24`

```python
def fret_chart(request):
    if is_ajax(request):
        forster_list = cache.get("forster_list")
        if not forster_list:
            job = cache.get("calc_fret_job")
            # ... celery job handling ...
            cache.set("forster_list", forster_list, 60 * 60 * 24)
```

- **What's cached:** Forster distance calculations
- **TTL:** 24 hours
- **Cache keys:** `"forster_list"`, `"calc_fret_job"`
- **Invalidation:** None

---

##### H. Forster List Calculation
**File:** `backend/proteins/util/helpers.py:321-414`

```python
def forster_list():
    # Try to get cached results first
    cache_key = "forster_list_results"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    # ... expensive calculation ...

    # Cache results for 6 hours
    cache.set(cache_key, result, 60 * 60 * 6)
    return result
```

- **What's cached:** Forster distance calculations for all protein pairs
- **TTL:** 6 hours
- **Cache key:** `"forster_list_results"`
- **Invalidation:** None
- **Note:** Includes memory optimization with batching and GC

---

##### I. MaxMind GeoIP Database
**File:** `backend/proteins/views/protein.py:125-188`

```python
def maxmind_db() -> str:
    cache_key = "maxmind_db_path"
    cached_path = cache.get(cache_key)
    if cached_path and os.path.exists(cached_path):
        return cached_path

    # Download and cache the database
    # ...
    cache.set(cache_key, tmp.name, 60 * 60 * 24)
    return tmp.name

def maxmind_reader() -> "maxminddb.Reader | None":
    cache_key = "maxmind_reader"
    reader = cache.get(cache_key)
    if reader is not None:
        return reader

    # ...
    cache.set(cache_key, reader, 60 * 60)
    return reader
```

- **What's cached:**
  - MaxMind DB file path (24 hours)
  - MaxMind DB reader object (1 hour)
- **Cache keys:** `"maxmind_db_path"`, `"maxmind_reader"`
- **Invalidation:** None (time-based only)
- **Issues:**
  - ⚠️ Caching file paths assumes file persistence

---

### 2. HTTP Caching Headers & View-Level Cache

#### 2.1 Page-Level Caching with @cache_page

##### A. Protein Detail Page
**File:** `backend/proteins/views/protein.py:220-223`

```python
class ProteinDetailView(DetailView):
    # Only enable caching in production (when DEBUG=False)
    if not settings.DEBUG:
        dispatch = method_decorator(cache_page(60 * 30))(dispatch)
        dispatch = method_decorator(vary_on_cookie)(dispatch)
```

- **What's cached:** Entire protein detail page HTML
- **TTL:** 30 minutes
- **Varies on:** Cookie (user-specific)
- **Invalidation:** ✅ Manual via `uncache_protein_page()` helper
- **Enabled:** Production only (`DEBUG=False`)

**Invalidation calls:**
- Line 419: After protein update
- Line 725: After adding reference
- Line 754: After adding excerpt
- Line 784: After reverting revision
- Line 814: After updating transitions
- Line 866: After updating bleach measurements
- Line 949: After flagging object

**Issues:**
- ⚠️ Manual invalidation is error-prone (must remember to call after every mutation)
- ⚠️ Easy to forget invalidation in new code paths

---

##### B. Spectra Image Generation
**File:** `backend/proteins/views/protein.py:531`

```python
@cache_page(60 * 120)
def spectra_image(request, slug, **kwargs):
```

- **What's cached:** Generated spectra images (PNG/SVG)
- **TTL:** 2 hours
- **Invalidation:** None
- **Issues:**
  - ⚠️ No invalidation when spectra change
  - ⚠️ Could serve stale images for up to 2 hours

---

##### C. Lineage Tree Data
**File:** `backend/proteins/views/ajax.py:240`

```python
@cache_page(60 * 5)
def get_lineage(request, slug=None, org=None):
```

- **What's cached:** Protein lineage tree JSON
- **TTL:** 5 minutes
- **Invalidation:** None
- **Issues:**
  - ⚠️ Short TTL may not be sufficient
  - ⚠️ No invalidation on lineage updates

---

##### D. Reference List
**File:** `backend/config/urls.py:121`

```python
path(
    "references/",
    cache_page(60 * 30)(ReferenceListView.as_view()),
    name="reference-list",
),
```

- **What's cached:** Reference list page
- **TTL:** 30 minutes
- **Invalidation:** None

---

##### E. REST API Endpoints
**File:** `backend/proteins/api/views.py`

**ProteinListAPIView2** (lines 62-64):
```python
@method_decorator(cache_page(60 * 10))
def dispatch(self, *args, **kwargs):
    return super().dispatch(*args, **kwargs)
```
- **TTL:** 10 minutes

**ProteinListAPIView** (lines 89-91):
```python
@method_decorator(cache_page(60 * 10))
def dispatch(self, *args, **kwargs):
    return super().dispatch(*args, **kwargs)
```
- **TTL:** 10 minutes

**ProteinTableAPIView** (lines 154-157):
```python
@method_decorator(cache_control(public=True, max_age=600))
@method_decorator(cache_page(60 * 10))
def dispatch(self, *args, **kwargs):
    return super().dispatch(*args, **kwargs)
```
- **TTL:** 10 minutes
- **HTTP headers:** ✅ `Cache-Control: public, max-age=600`
- **Issues:**
  - ⚠️ No Vary headers (all users see same cached data)
  - ⚠️ No ETags

---

#### 2.2 Manual Cache Invalidation Helper

**File:** `backend/fpbase/util.py:1-43`

```python
def get_view_cache_key(view_name, args=None, namespace=None, key_prefix=None, request=None):
    """Get cache key for view-level cache"""
    # ...
    req.path = reverse(view_name, args=args)
    return get_cache_key(req, key_prefix=key_prefix)

def clear_view_cache(*args, **kwargs):
    cache_key = get_view_cache_key(*args, **kwargs)
    key_deleted = False
    if cache_key:
        key_deleted = cache.delete(cache_key)
    return key_deleted

def uncache_protein_page(slug, request):
    clear_view_cache("proteins:protein-detail", args=[slug], request=request)
```

- **Purpose:** Manually invalidate `@cache_page` decorated views
- **Used for:** Protein detail pages after updates
- **Issues:**
  - 🔴 **CRITICAL:** Manual invalidation is error-prone
  - ⚠️ Must remember to call after all mutations
  - ⚠️ Fragile - breaks if cache key generation changes

---

### 3. Template Caching

**Status:** ❌ **NOT FOUND**
- No `{% cache %}` template tags found in the codebase
- Opportunity: Could add fragment caching for expensive template sections

---

### 4. Database Query Caching

#### 4.1 Python cached_property

**Spectrum Model** (`backend/proteins/models/spectrum.py:598-606`):
```python
@cached_property
def x(self):
    """Extract x values from data. Cached to avoid repeated allocations."""
    return [i[0] for i in self.data]

@cached_property
def y(self):
    """Extract y values from data. Cached to avoid repeated allocations."""
    return [i[1] for i in self.data]
```

- **What's cached:** Extracted x/y arrays from spectrum data
- **Scope:** Per-instance (in-memory)
- **Invalidation:** ✅ Manually cleared in `change_x()`/`change_y()` methods
- **Issues:** None (appropriate use of cached_property)

**Microscope Model** (`backend/proteins/models/microscope.py:90-158`):
```python
@cached_property
def has_inverted_bs(self): ...

@cached_property
def has_reflective_emfilters(self): ...

@cached_property
def lights(self): ...

@cached_property
def cameras(self): ...

@cached_property
def lasers(self): ...

@cached_property
def ex_filters(self): ...

@cached_property
def em_filters(self): ...

@cached_property
def bs_filters(self): ...

@cached_property
def spectra(self): ...
```

- **What's cached:** Various microscope-related querysets and computations
- **Scope:** Per-instance (in-memory)
- **Invalidation:** None
- **Issues:**
  - ⚠️ Stale data if DB changes after object is loaded
  - ⚠️ Should document as "instance-scoped only"

---

### 5. GraphQL Caching

**File:** `backend/proteins/schema/query.py`

**Patterns:**
- Uses `get_cached_spectrum()` (see section 1.2.B - caches ORM objects ❌)
- Uses spectra sluglist caching (see section 1.2.B)
- No GraphQL-specific caching layer (relies on Django cache)
- No response-level HTTP caching headers

---

### 6. REST API Caching

**File:** `backend/proteins/api/views.py`

**Patterns:**
- View-level caching with `@cache_page` (see section 2.1.E)
- HTTP `Cache-Control` headers on `ProteinTableAPIView` only
- No DRF-specific caching (throttling only)
- Missing: Vary headers, ETags, Last-Modified

---

### 7. Proxy/CDN Caching

#### 7.1 AWS S3 Static Files
**File:** `backend/config/settings/production.py:93-97`

```python
_AWS_EXPIRY = 60 * 60 * 24 * 7
AWS_S3_OBJECT_PARAMETERS = {
    "CacheControl": f"max-age={_AWS_EXPIRY}, s-maxage={_AWS_EXPIRY}, must-revalidate",
}
```

- **What's cached:** Uploaded media files on S3
- **TTL:** 7 days (`max-age=604800, s-maxage=604800`)
- **HTTP header:** ✅ `Cache-Control: max-age=604800, s-maxage=604800, must-revalidate`
- **CDN behavior:** Both browser and CDN cache for 7 days
- **Issues:** None (appropriate for user-uploaded media)

---

#### 7.2 WhiteNoise Static Files
**File:** `backend/config/settings/production.py:115-122`

```python
WHITENOISE_MAX_AGE = 600

def WHITENOISE_IMMUTABLE_FILE_TEST(path, url):
    # Match vite (rollup)-generated hashes
    return re.match(r"^.+[.-][0-9a-zA-Z_-]{8,12}\..+$", url)
```

- **What's cached:** Static assets (CSS, JS, images)
- **TTL:**
  - Default: 10 minutes (`600`)
  - Immutable (hashed) files: Far-future expiry
- **HTTP headers:** ✅ Set by WhiteNoise based on file hash
- **Issues:** None (WhiteNoise handles this correctly)

---

#### 7.3 Cloudflare CDN

**Status:** ⚠️ **INFERRED (NOT EXPLICITLY CONFIGURED IN CODEBASE)**

- Found 34 files mentioning "Cloudflare" (mostly frontend JS, templates)
- No explicit Cloudflare configuration in Django settings
- Likely configured at DNS/proxy level (not in codebase)
- **No page rules, cache rules, or Cloudflare-specific headers found**
- **Cannot programmatically purge cache** (no API integration)

---

## Cache Key Inventory

| Cache Key | TTL | Purpose | Invalidation | File |
|-----------|-----|---------|--------------|------|
| `gbseqs` | 24h | GenBank sequences | None | entrez.py:144 |
| `_spectrum_{id}` | 24h | Individual spectra (ORM objects ❌) | None | query.py:30 |
| `_spectra_sluglist_{filters}` | 1h | Filtered spectra list | None | query.py:116 |
| `spectra_sluglist` | 1h | All spectra info | ✅ Manual on save | spectrum.py:79, 424 |
| `optical_config_list` | 1h | Optical configs | ✅ Manual on save | microscope.py:198, 233 |
| `slug_dict` | **60s** 🔴 | Protein name→slug map | None | helpers.py:42 |
| `ga_popular_proteins` | 24h | Google Analytics data | None | ga.py:54 |
| `forster_list` | 24h | FRET calculations | None | fret.py:24 |
| `forster_list_results` | 6h | Forster distances | None | helpers.py:413 |
| `calc_fret_job` | N/A | Celery job ID | Manual | fret.py:29 |
| `maxmind_db_path` | 24h | GeoIP DB file path | None | protein.py:160 |
| `maxmind_reader` | 1h | GeoIP DB reader | None | protein.py:183 |

---

## Critical Issues Summary

### 🔴 Critical Issues (Fix First)

1. **Caching ORM objects** (`get_cached_spectrum()`)
   - Can cause stale relationships, Django version incompatibility
   - **Fix:** Cache serialized JSON instead, or remove caching

2. **60-second slug_dict TTL**
   - Causes frequent cache misses on every page
   - **Fix:** Increase to 1 hour + signal-based invalidation

3. **No ETags on large payloads** (Issue #361)
   - 140 KB re-downloaded every 10 minutes even when unchanged
   - **Fix:** Add ETag support to spectra endpoints

4. **Manual cache invalidation is error-prone**
   - Easy to forget `uncache_protein_page()` calls
   - **Fix:** Signal-based automatic invalidation

### 🟡 Medium Issues

5. **Unbounded cache growth** (`gbseqs`)
   - Dictionary grows forever, never clears old entries
   - **Fix:** Individual cache keys per GenBank ID

6. **No Vary headers on API endpoints**
   - CDN/browsers can serve wrong cached data
   - **Fix:** Add `Vary: Accept-Encoding, Cookie` where appropriate

7. **Caching file paths** (MaxMind)
   - Assumes temporary files persist
   - **Fix:** Use named temp files or store in media directory

8. **Arbitrary TTLs everywhere**
   - Magic numbers (60, 300, 3600, 86400) with no strategy
   - **Fix:** Centralized TTL constants + semantic naming

### 🟢 Minor Issues

9. **No template fragment caching**
   - Opportunity for performance gains
   - **Consider:** Add `{% cache %}` for expensive sections

10. **cached_property not invalidated**
    - Stale data if DB changes after object load
    - **Fix:** Document as instance-scoped, or use django-cache-memoize

11. **No CDN purge integration**
    - Can't programmatically clear Cloudflare cache
    - **Defer:** Would need API integration (not in Django codebase)

---

## Implementation Plan

### Phase 1: Infrastructure & Version Tracking System

#### 1.1 Create Cache Versioning Infrastructure
**New file:** `backend/fpbase/cache.py`

```python
"""
Cache versioning and invalidation infrastructure.

Provides version-based cache keys for automatic invalidation and ETag support.
"""
from django.core.cache import cache
from typing import Optional


class CacheVersion:
    """Manages version numbers for cache keys to enable invalidation."""

    VERSION_PREFIX = "cache_version:"

    @classmethod
    def get_version(cls, key: str) -> int:
        """Get current version number for a cache key."""
        version_key = f"{cls.VERSION_PREFIX}{key}"
        version = cache.get(version_key)
        if version is None:
            version = 1
            cache.set(version_key, version, timeout=None)  # Never expire
        return version

    @classmethod
    def increment_version(cls, key: str) -> int:
        """Increment version number (invalidates all caches using this key)."""
        version_key = f"{cls.VERSION_PREFIX}{key}"
        try:
            new_version = cache.incr(version_key)
        except ValueError:
            # Key doesn't exist, initialize it
            new_version = 1
            cache.set(version_key, new_version, timeout=None)
        return new_version

    @classmethod
    def generate_etag(cls, key: str) -> str:
        """Generate a strong ETag from version number."""
        version = cls.get_version(key)
        return f'"{version}"'

    @classmethod
    def versioned_key(cls, base_key: str, resource_type: str) -> str:
        """Create a cache key that includes the version number."""
        version = cls.get_version(resource_type)
        return f"{base_key}:v{version}"


# Resource type constants
RESOURCE_SPECTRA = "spectra"
RESOURCE_PROTEINS = "proteins"
RESOURCE_OPTICAL_CONFIGS = "optical_configs"
RESOURCE_REFERENCES = "references"
```

**Implementation tasks:**
- [ ] Create `backend/fpbase/cache.py` with `CacheVersion` class
- [ ] Add thread-safe atomic increment operations
- [ ] Support multiple resource types (spectra, proteins, optical_configs, etc.)
- [ ] Add utility methods for versioned cache keys
- [ ] Write unit tests for version tracking

---

#### 1.2 Create Signal-Based Cache Invalidation
**New file:** `backend/fpbase/signals.py`

```python
"""
Django signals for automatic cache invalidation.

Connects to model save/delete signals to invalidate caches automatically.
"""
from django.db.models.signals import post_save, post_delete, m2m_changed
from django.dispatch import receiver
from django.core.cache import cache

from proteins.models import Spectrum, Protein, OpticalConfig, Reference
from .cache import CacheVersion, RESOURCE_SPECTRA, RESOURCE_PROTEINS, RESOURCE_OPTICAL_CONFIGS


@receiver([post_save, post_delete], sender=Spectrum)
def invalidate_spectrum_cache(sender, instance, **kwargs):
    """Invalidate spectrum caches when a spectrum is saved or deleted."""
    # Increment version (invalidates all spectrum-related caches)
    CacheVersion.increment_version(RESOURCE_SPECTRA)

    # Clear individual spectrum cache
    cache.delete(f"_spectrum_{instance.id}")


@receiver([post_save, post_delete], sender=Protein)
def invalidate_protein_cache(sender, instance, **kwargs):
    """Invalidate protein caches when a protein is saved or deleted."""
    CacheVersion.increment_version(RESOURCE_PROTEINS)

    # Clear slug dictionary (protein names changed)
    cache.delete("slug_dict")

    # Clear forster calculations (protein properties changed)
    cache.delete("forster_list")
    cache.delete("forster_list_results")


@receiver([post_save, post_delete], sender=OpticalConfig)
def invalidate_optical_config_cache(sender, instance, **kwargs):
    """Invalidate optical config caches."""
    CacheVersion.increment_version(RESOURCE_OPTICAL_CONFIGS)


@receiver([post_save, post_delete], sender=Reference)
def invalidate_reference_cache(sender, instance, **kwargs):
    """Invalidate reference caches."""
    # If reference is linked to a protein, invalidate that protein's page
    if hasattr(instance, 'protein') and instance.protein:
        CacheVersion.increment_version(RESOURCE_PROTEINS)
```

**Implementation tasks:**
- [ ] Create `backend/fpbase/signals.py`
- [ ] Connect signals for all relevant models
- [ ] Replace all manual `cache.delete()` calls with version increments
- [ ] Add `apps.py` to register signals on app startup
- [ ] Test signal firing on save/delete/m2m_changed

---

#### 1.3 Create HTTP Caching Utilities
**New file:** `backend/fpbase/http_cache.py`

```python
"""
HTTP caching utilities for ETags and Cache-Control headers.
"""
from functools import wraps
from django.http import HttpResponse
from django.utils.cache import patch_cache_control, patch_vary_headers
from .cache import CacheVersion


def conditional_get(resource_type):
    """
    Decorator to add ETag support and handle conditional GET requests.

    Returns 304 Not Modified if If-None-Match header matches current ETag.

    Usage:
        @conditional_get(RESOURCE_SPECTRA)
        def my_view(request):
            return JsonResponse(data)
    """
    def decorator(view_func):
        @wraps(view_func)
        def wrapper(request, *args, **kwargs):
            # Generate ETag from resource version
            etag = CacheVersion.generate_etag(resource_type)

            # Check If-None-Match header
            if_none_match = request.META.get('HTTP_IF_NONE_MATCH')
            if if_none_match == etag:
                # Client has current version, return 304
                response = HttpResponse(status=304)
                response['ETag'] = etag
                return response

            # Client needs new data, call view
            response = view_func(request, *args, **kwargs)

            # Add ETag to response
            if response.status_code == 200:
                response['ETag'] = etag

            return response
        return wrapper
    return decorator


def cache_control_public(max_age=600):
    """
    Decorator to set Cache-Control: public, max-age=X header.

    Use for public data that can be cached by browsers and CDNs.
    """
    def decorator(view_func):
        @wraps(view_func)
        def wrapper(request, *args, **kwargs):
            response = view_func(request, *args, **kwargs)
            patch_cache_control(response, public=True, max_age=max_age)
            return response
        return wrapper
    return decorator


def cache_control_private(max_age=600):
    """
    Decorator to set Cache-Control: private, max-age=X header.

    Use for user-specific data that should only be cached by browsers.
    """
    def decorator(view_func):
        @wraps(view_func)
        def wrapper(request, *args, **kwargs):
            response = view_func(request, *args, **kwargs)
            patch_cache_control(response, private=True, max_age=max_age)
            return response
        return wrapper
    return decorator


def vary_on_headers(*headers):
    """
    Decorator to add Vary headers.

    Use to tell caches that the response varies based on request headers.

    Usage:
        @vary_on_headers('Accept-Encoding', 'Cookie')
        def my_view(request):
            return JsonResponse(data)
    """
    def decorator(view_func):
        @wraps(view_func)
        def wrapper(request, *args, **kwargs):
            response = view_func(request, *args, **kwargs)
            patch_vary_headers(response, headers)
            return response
        return wrapper
    return decorator
```

**Implementation tasks:**
- [ ] Create `backend/fpbase/http_cache.py`
- [ ] Implement `@conditional_get` decorator for ETag support
- [ ] Implement `@cache_control_public` and `@cache_control_private` decorators
- [ ] Implement `@vary_on_headers` decorator
- [ ] Add integration tests with actual HTTP requests

---

### Phase 2: Fix Spectra Caching (Issue #361)

#### 2.1 Add ETag Support to REST API
**File:** `backend/proteins/api/views.py`

**Before:**
```python
@require_http_methods(["GET", "HEAD"])
def spectraslugs(request):
    spectrainfo = get_cached_spectra_info()
    response = HttpResponse(spectrainfo, content_type="application/json")
    response["Cache-Control"] = "public, max-age=600"
    return response
```

**After:**
```python
from fpbase.http_cache import conditional_get, cache_control_public, vary_on_headers
from fpbase.cache import RESOURCE_SPECTRA

@require_http_methods(["GET", "HEAD"])
@conditional_get(RESOURCE_SPECTRA)
@cache_control_public(max_age=600)
@vary_on_headers('Accept-Encoding')
def spectraslugs(request):
    spectrainfo = get_cached_spectra_info()
    response = HttpResponse(spectrainfo, content_type="application/json")
    return response
```

**Expected behavior:**
```http
# First request
GET /api/proteins/spectraslugs/
→ 200 OK
  ETag: "123"
  Cache-Control: public, max-age=600
  Vary: Accept-Encoding
  [140 KB payload]

# After 10+ minutes (cache expired)
GET /api/proteins/spectraslugs/
If-None-Match: "123"
→ 304 Not Modified
  ETag: "123"
  [NO PAYLOAD - just headers, < 1 KB!]

# After spectrum is saved (version incremented)
GET /api/proteins/spectraslugs/
If-None-Match: "123"
→ 200 OK
  ETag: "124"  ← New version!
  [140 KB new payload]
```

**Implementation tasks:**
- [ ] Update `spectraslugs()` view with decorators
- [ ] Add `Vary: Accept-Encoding` header
- [ ] Test 304 responses with If-None-Match
- [ ] Measure bandwidth savings in production

---

#### 2.2 Add ETag Support to GraphQL
**File:** `backend/fpbase/views.py`

```python
from fpbase.cache import CacheVersion, RESOURCE_SPECTRA

class RateLimitedGraphQLView(GraphQLView):
    def dispatch(self, request, *args, **kwargs):
        # Check if this is a SpectraList query
        is_spectra_query = False
        if request.body and b'SpectraList' in request.body:
            is_spectra_query = True

            # Generate ETag from spectra version
            etag = CacheVersion.generate_etag(RESOURCE_SPECTRA)

            # Check If-None-Match header
            if_none_match = request.META.get('HTTP_IF_NONE_MATCH')
            if if_none_match == etag:
                # Client has current version
                response = HttpResponse(status=304)
                response['ETag'] = etag
                response['Cache-Control'] = 'public, max-age=600'
                return response

        # Process request normally
        response = super().dispatch(request, *args, **kwargs)

        # Add ETag to successful SpectraList responses
        if is_spectra_query and response.status_code == 200:
            response['ETag'] = etag
            response['Cache-Control'] = 'public, max-age=600'

        return response
```

**Implementation tasks:**
- [ ] Add ETag detection for SpectraList queries
- [ ] Return 304 for conditional requests
- [ ] Test with GraphQL client
- [ ] Ensure backward compatibility

---

#### 2.3 Update Spectrum Model Signal
**File:** `backend/proteins/models/spectrum.py`

**Before:**
```python
def save(self, *args, **kwargs):
    cache.delete(SPECTRA_CACHE_KEY)
    super().save(*args, **kwargs)
```

**After:**
```python
def save(self, *args, **kwargs):
    # Signal handler will automatically increment version
    # See backend/fpbase/signals.py
    super().save(*args, **kwargs)
```

**Implementation tasks:**
- [ ] Remove manual `cache.delete()` from `save()` method
- [ ] Verify signal-based invalidation works
- [ ] Test that ETag increments on save
- [ ] Cleanup unused `SPECTRA_CACHE_KEY` constant

---

### Phase 3: Fix Critical Caching Issues

#### 3.1 Stop Caching ORM Objects
**File:** `backend/proteins/schema/query.py`

**Problem:** Caching ORM objects is problematic:
- Relationships can become stale
- Django version upgrades can break pickled objects
- Wastes memory (storing entire objects)

**Option A: Cache serialized JSON**
```python
def get_cached_spectrum_data(id, timeout=60 * 60 * 24):
    """Cache spectrum data as JSON, not ORM objects."""
    key = f"_spectrum_json_{id}"
    data = cache.get(key)
    if not data:
        try:
            spectrum = (
                models.Spectrum.objects.filter(id=id)
                .select_related(...)
                .get()
            )
            # Serialize to JSON-compatible dict
            data = {
                'id': spectrum.id,
                'name': spectrum.name,
                'owner_id': spectrum.owner_id if spectrum.owner else None,
                # ... other fields
            }
            cache.set(key, data, timeout)
        except models.Spectrum.DoesNotExist:
            return None
    return data
```

**Option B: Remove caching entirely**
```python
def get_spectrum(id):
    """No caching - PostgreSQL is fast enough for single-object lookups."""
    try:
        return (
            models.Spectrum.objects.filter(id=id)
            .select_related(...)
            .get()
        )
    except models.Spectrum.DoesNotExist:
        return None
```

**Recommendation:** Try Option B first (measure performance), fall back to Option A if needed.

**Implementation tasks:**
- [ ] Benchmark current performance with caching
- [ ] Implement Option B (remove caching)
- [ ] Measure performance without caching
- [ ] If needed, implement Option A (JSON caching)
- [ ] Update all call sites

---

#### 3.2 Fix Unbounded Cache Growth
**File:** `backend/proteins/extrest/entrez.py`

**Before:**
```python
def get_cached_gbseqs(gbids, max_age=60 * 60 * 24):
    gbseqs = cache.get("gbseqs", {})  # Unbounded dict
    now = time.time()
    tofetch = [id for id in gbids if id not in gbseqs or gbseqs[id][1] - now >= max_age]
    gbseqs.update({k: (v, now) for k, v in _fetch_gb_seqs(tofetch).items()})
    cache.set("gbseqs", gbseqs, 60 * 60 * 24)  # Store entire dict
    return gbseqs
```

**After:**
```python
def get_cached_gbseqs(gbids, max_age=60 * 60 * 24):
    """Fetch GenBank sequences, using individual cache keys per ID."""
    results = {}
    tofetch = []

    # Check cache for each ID individually
    for gbid in gbids:
        cache_key = f"gbseq:{gbid}"
        cached = cache.get(cache_key)
        if cached is not None:
            results[gbid] = cached
        else:
            tofetch.append(gbid)

    # Fetch missing sequences
    if tofetch:
        fetched = _fetch_gb_seqs(tofetch)
        for gbid, seq in fetched.items():
            cache.set(f"gbseq:{gbid}", seq, max_age)
            results[gbid] = seq

    return results
```

**Benefits:**
- Each GenBank ID has its own TTL
- No unbounded dictionary growth
- Old entries expire naturally
- No manual timestamp tracking

**Implementation tasks:**
- [ ] Rewrite `get_cached_gbseqs()` to use individual keys
- [ ] Remove manual timestamp tracking
- [ ] Test with multiple GenBank IDs
- [ ] Verify old cache entries expire correctly

---

#### 3.3 Fix slug_dict Short TTL
**File:** `backend/proteins/util/helpers.py`

**Before:**
```python
slug_dict = cache.get_or_set("slug_dict", create_slug_dict, 60)  # 60 seconds!
```

**After:**
```python
from fpbase.cache import CacheVersion, RESOURCE_PROTEINS

def get_slug_dict():
    """Get protein slug dictionary with version-based invalidation."""
    version = CacheVersion.get_version(RESOURCE_PROTEINS)
    cache_key = f"slug_dict:v{version}"

    slug_dict = cache.get(cache_key)
    if slug_dict is None:
        slug_dict = create_slug_dict()
        cache.set(cache_key, slug_dict, 60 * 60)  # 1 hour

    return slug_dict
```

**Changes:**
- TTL increased from 60 seconds to 1 hour
- Version-based cache key (auto-invalidates when proteins change)
- Signal handler increments version on Protein save/delete

**Implementation tasks:**
- [ ] Create `get_slug_dict()` helper function
- [ ] Update all call sites to use new function
- [ ] Add signal handler for Protein model (already in Phase 1.2)
- [ ] Test invalidation when protein is saved
- [ ] Measure cache hit rate improvement

---

#### 3.4 Replace Manual Page Cache Invalidation
**File:** `backend/proteins/views/protein.py`

**Before:**
```python
# In 7+ different places:
uncache_protein_page(slug, request)
```

**After:**
```python
# Option A: Version-based cache keys (preferred)
class ProteinDetailView(DetailView):
    if not settings.DEBUG:
        @method_decorator(cache_page_versioned(60 * 30, RESOURCE_PROTEINS))
        @method_decorator(vary_on_cookie)
        def dispatch(self, *args, **kwargs):
            return super().dispatch(*args, **kwargs)

# Option B: Remove page caching entirely, rely on HTTP caching
class ProteinDetailView(DetailView):
    # No @cache_page decorator
    # Just set Cache-Control and ETag headers in response
```

**Recommendation:** Option B (remove page caching, use HTTP caching only)

**Implementation tasks:**
- [ ] Remove `@cache_page` decorator from ProteinDetailView
- [ ] Add `@conditional_get(RESOURCE_PROTEINS)` decorator
- [ ] Add proper Cache-Control headers
- [ ] Remove all `uncache_protein_page()` calls
- [ ] Delete `backend/fpbase/util.py` (no longer needed)
- [ ] Test that pages update immediately after protein changes

---

### Phase 4: Add Proper HTTP Headers

#### 4.1 Add Vary Headers to API Endpoints
**File:** `backend/proteins/api/views.py`

**Add to all public API views:**
```python
from fpbase.http_cache import vary_on_headers

@method_decorator(vary_on_headers('Accept-Encoding'))
def dispatch(self, *args, **kwargs):
    return super().dispatch(*args, **kwargs)
```

**Add to user-specific views:**
```python
@method_decorator(vary_on_headers('Accept-Encoding', 'Cookie'))
def dispatch(self, *args, **kwargs):
    return super().dispatch(*args, **kwargs)
```

**Why this matters:**
- `Vary: Accept-Encoding` - CDN caches gzipped and non-gzipped versions separately
- `Vary: Cookie` - CDN doesn't serve authenticated user's data to anonymous users

**Implementation tasks:**
- [ ] Add `Vary: Accept-Encoding` to all API views
- [ ] Add `Vary: Cookie` to user-specific endpoints
- [ ] Test CDN respects Vary headers
- [ ] Document which endpoints vary on what

---

#### 4.2 Standardize Cache-Control Headers
**Create:** `backend/fpbase/cache_config.py`

```python
"""
Centralized cache TTL configuration.

Use these semantic constants instead of magic numbers.
"""

# HTTP cache TTLs (browser/CDN)
TTL_5_MIN = 60 * 5      # Frequently changing data
TTL_10_MIN = 60 * 10    # Default for most API endpoints
TTL_1_HOUR = 60 * 60    # Stable data
TTL_1_DAY = 60 * 60 * 24    # Very stable data
TTL_1_WEEK = 60 * 60 * 24 * 7   # Media files

# Server-side cache TTLs (Redis)
CACHE_TTL_SHORT = 60 * 5    # 5 minutes
CACHE_TTL_MEDIUM = 60 * 60  # 1 hour
CACHE_TTL_LONG = 60 * 60 * 6    # 6 hours
CACHE_TTL_VERY_LONG = 60 * 60 * 24  # 24 hours

# Cache key prefixes
CACHE_KEY_SPECTRUM = "spectrum"
CACHE_KEY_PROTEIN = "protein"
CACHE_KEY_OPTICAL_CONFIG = "optical_config"
CACHE_KEY_GENBANK = "gbseq"
CACHE_KEY_FORSTER = "forster"
CACHE_KEY_SLUG_DICT = "slug_dict"
CACHE_KEY_GA_POPULAR = "ga_popular"
```

**Update all cache calls to use constants:**
```python
# Before
cache.set("slug_dict", data, 60)

# After
from fpbase.cache_config import CACHE_TTL_MEDIUM, CACHE_KEY_SLUG_DICT
cache.set(CACHE_KEY_SLUG_DICT, data, CACHE_TTL_MEDIUM)
```

**Implementation tasks:**
- [ ] Create `backend/fpbase/cache_config.py`
- [ ] Define semantic TTL constants
- [ ] Replace all magic numbers with constants
- [ ] Document when to use each TTL
- [ ] Add to code review checklist

---

#### 4.3 Add ETag Support to Common Views

**Target views:**
- ✅ Spectra sluglist (Phase 2)
- Protein detail page
- Spectra image generation
- Lineage tree JSON
- Reference list
- Optical config list

**Pattern (repeat for each view):**
```python
from fpbase.http_cache import conditional_get, cache_control_public, vary_on_headers
from fpbase.cache import RESOURCE_PROTEINS

@conditional_get(RESOURCE_PROTEINS)
@cache_control_public(max_age=600)
@vary_on_headers('Accept-Encoding')
def my_view(request):
    # ... generate response ...
    return response
```

**Implementation tasks:**
- [ ] Add ETags to spectra image view
- [ ] Add ETags to lineage tree view
- [ ] Add ETags to reference list view
- [ ] Add ETags to optical config view
- [ ] Measure bandwidth reduction for each
- [ ] Document expected cache hit rates

---

### Phase 5: Cleanup & Documentation

#### 5.1 Remove Dead Code

**Files to delete:**
- [ ] `backend/fpbase/util.py` (manual cache invalidation helpers)

**Code to remove:**
- [ ] All `uncache_protein_page()` calls
- [ ] All manual `cache.delete()` calls (replaced by signals)
- [ ] Unused cache key constants

**Implementation tasks:**
- [ ] Search for `uncache_protein_page` and remove all calls
- [ ] Search for `cache.delete` and replace with version increments
- [ ] Run tests to ensure nothing breaks
- [ ] Clean up imports

---

#### 5.2 Create Cache Configuration Documentation

**Update:** `README.md` (add new section)

```markdown
## Caching Strategy

FPbase uses a multi-layered caching strategy:

### 1. HTTP Caching (ETags + Cache-Control)
- **What:** Browser and CDN caching via HTTP headers
- **TTL:** 10 minutes (default), up to 1 hour for stable data
- **Invalidation:** ETags enable conditional requests (304 Not Modified)
- **Use for:** API responses, public pages, images

### 2. Server-Side Caching (Redis)
- **What:** Django cache framework with Redis backend
- **TTL:** 5 minutes to 24 hours depending on data stability
- **Invalidation:** Automatic via Django signals (version-based)
- **Use for:** Expensive computations, external API calls

### 3. Application Caching (cached_property)
- **What:** In-memory per-instance caching
- **TTL:** Lifetime of Python object
- **Invalidation:** None (instance-scoped only)
- **Use for:** Expensive property calculations

### How to Add Caching

#### For API Endpoints (Recommended Pattern)
\`\`\`python
from fpbase.http_cache import conditional_get, cache_control_public, vary_on_headers
from fpbase.cache import RESOURCE_SPECTRA

@conditional_get(RESOURCE_SPECTRA)  # Add ETag support
@cache_control_public(max_age=600)   # 10 minute browser cache
@vary_on_headers('Accept-Encoding')  # CDN varies on encoding
def my_api_view(request):
    data = expensive_calculation()
    return JsonResponse(data)
\`\`\`

#### For Expensive Computations
\`\`\`python
from django.core.cache import cache
from fpbase.cache_config import CACHE_TTL_MEDIUM

def expensive_calculation():
    cache_key = "my_calculation"
    result = cache.get(cache_key)
    if result is None:
        result = do_expensive_work()
        cache.set(cache_key, result, CACHE_TTL_MEDIUM)
    return result
\`\`\`

### Cache Invalidation

**Automatic (preferred):**
Cache invalidation happens automatically via Django signals. When you save/delete a model, the version number increments and all related caches are invalidated.

**Manual (rare):**
Only needed for external data sources:
\`\`\`python
from fpbase.cache import CacheVersion, RESOURCE_SPECTRA
CacheVersion.increment_version(RESOURCE_SPECTRA)
\`\`\`

### Cache TTL Guidelines

| Data Type | TTL | Rationale |
|-----------|-----|-----------|
| Frequently changing (e.g., GA stats) | 5 minutes | Fresh data important |
| API endpoints (default) | 10 minutes | Balance freshness & performance |
| Stable data (e.g., protein properties) | 1 hour | Changes infrequent |
| External APIs (e.g., GenBank) | 24 hours | Avoid rate limits |
| Static media | 7 days | Content-addressed (hashed URLs) |

See `backend/fpbase/cache_config.py` for TTL constants.
```

---

#### 5.3 Add Tests
**New file:** `backend/tests/test_cache.py`

```python
"""
Tests for caching infrastructure.
"""
import pytest
from django.core.cache import cache
from django.test import RequestFactory
from fpbase.cache import CacheVersion, RESOURCE_SPECTRA
from proteins.models import Spectrum


class TestCacheVersion:
    def test_get_version_initializes(self):
        """Version starts at 1 if not set."""
        cache.clear()
        version = CacheVersion.get_version("test_resource")
        assert version == 1

    def test_increment_version(self):
        """Incrementing version increases by 1."""
        cache.clear()
        v1 = CacheVersion.get_version("test_resource")
        v2 = CacheVersion.increment_version("test_resource")
        assert v2 == v1 + 1

    def test_generate_etag(self):
        """ETag format is correct."""
        cache.clear()
        etag = CacheVersion.generate_etag("test_resource")
        assert etag == '"1"'  # Strong ETag

    def test_versioned_key(self):
        """Versioned cache keys include version number."""
        cache.clear()
        key = CacheVersion.versioned_key("data", "test_resource")
        assert key == "data:v1"

        CacheVersion.increment_version("test_resource")
        key2 = CacheVersion.versioned_key("data", "test_resource")
        assert key2 == "data:v2"


class TestSpectrumSignals:
    def test_spectrum_save_increments_version(self, db):
        """Saving a spectrum increments the spectra version."""
        cache.clear()
        initial_version = CacheVersion.get_version(RESOURCE_SPECTRA)

        spectrum = Spectrum.objects.create(...)

        new_version = CacheVersion.get_version(RESOURCE_SPECTRA)
        assert new_version == initial_version + 1

    def test_spectrum_delete_increments_version(self, db):
        """Deleting a spectrum increments the spectra version."""
        spectrum = Spectrum.objects.create(...)
        cache.clear()
        initial_version = CacheVersion.get_version(RESOURCE_SPECTRA)

        spectrum.delete()

        new_version = CacheVersion.get_version(RESOURCE_SPECTRA)
        assert new_version == initial_version + 1


class TestConditionalGet:
    def test_etag_in_response(self, client):
        """API responses include ETag header."""
        response = client.get('/api/proteins/spectraslugs/')
        assert 'ETag' in response
        assert response['ETag'].startswith('"')

    def test_304_on_matching_etag(self, client):
        """Returns 304 when If-None-Match matches ETag."""
        # First request
        response1 = client.get('/api/proteins/spectraslugs/')
        etag = response1['ETag']

        # Second request with If-None-Match
        response2 = client.get(
            '/api/proteins/spectraslugs/',
            HTTP_IF_NONE_MATCH=etag
        )
        assert response2.status_code == 304
        assert len(response2.content) == 0  # No body

    def test_200_on_different_etag(self, client, db):
        """Returns 200 when ETag changes."""
        # First request
        response1 = client.get('/api/proteins/spectraslugs/')
        etag1 = response1['ETag']

        # Change data (increments version)
        Spectrum.objects.create(...)

        # Second request with old ETag
        response2 = client.get(
            '/api/proteins/spectraslugs/',
            HTTP_IF_NONE_MATCH=etag1
        )
        assert response2.status_code == 200
        assert response2['ETag'] != etag1
        assert len(response2.content) > 0
```

**Implementation tasks:**
- [ ] Create `backend/tests/test_cache.py`
- [ ] Test version tracking
- [ ] Test signal-based invalidation
- [ ] Test ETag generation
- [ ] Test 304 responses
- [ ] Test Cache-Control headers
- [ ] Achieve >90% coverage for caching code

---

#### 5.4 Update Documentation

**Files to update:**
- [ ] `README.md` - Add caching strategy section (see 5.2)
- [ ] `CONTRIBUTING.md` - Add guidelines for when/how to add caching
- [ ] Code comments - Document cache invalidation flow

**Create new docs:**
- [ ] `docs/caching.md` - Detailed technical documentation
- [ ] `docs/performance.md` - Performance benchmarks before/after

---

## Implementation Order (Recommended)

### Sprint 1: Foundation (1-2 days)
1. ✅ **Phase 1.1** - Create cache versioning infrastructure
2. ✅ **Phase 1.2** - Create signal-based invalidation
3. ✅ **Phase 1.3** - Create HTTP caching utilities
4. ✅ **Phase 5.3** - Add tests for infrastructure

### Sprint 2: Quick Wins (1-2 days)
5. ✅ **Phase 3.3** - Fix slug_dict short TTL
6. ✅ **Phase 4.2** - Standardize Cache-Control headers
7. ✅ **Phase 3.1** - Stop caching ORM objects

### Sprint 3: Issue #361 (1 day)
8. ✅ **Phase 2.1** - Add ETag to REST API
9. ✅ **Phase 2.2** - Add ETag to GraphQL
10. ✅ **Phase 2.3** - Update Spectrum model signal

### Sprint 4: Polish & Extend (2-3 days)
11. ✅ **Phase 4.3** - Add ETags to other views
12. ✅ **Phase 4.1** - Add Vary headers
13. ✅ **Phase 3.4** - Replace manual page cache invalidation
14. ✅ **Phase 3.2** - Fix GenBank cache

### Sprint 5: Cleanup (1 day)
15. ✅ **Phase 5.1** - Remove dead code
16. ✅ **Phase 5.2** - Update documentation
17. ✅ **Phase 5.4** - Final documentation pass

**Total estimated time:** 6-9 days (1-2 weeks)

---

## Design Decisions to Confirm

### 1. ETag Format
**Question:** Strong ETags (`"v123"`) or weak (`W/"v123"`)?

**Recommendation:** **Strong ETags**
- Data is deterministic (same version = identical response)
- Enables byte-range requests (future optimization)
- More cacheable by intermediaries

---

### 2. Version Storage
**Question:** Store versions in Redis (existing cache) or database table?

**Recommendation:** **Redis**
- Already configured and running
- Fast atomic increments (`cache.incr()`)
- No DB migrations needed
- Versions don't need persistence (can rebuild if Redis flushes)

**Alternative:** If Redis is flushed, versions reset to 1 (clients re-download once, then ETags work again)

---

### 3. GraphQL SpectraList Caching
**Question:** Add ETags to GraphQL, migrate to REST, or both?

**Recommendation:** **Both (Option C from issue #361)**
- Add ETags to GraphQL (backward compatible)
- Keep REST endpoint with ETags
- Let clients choose (GraphQL for flexibility, REST for simplicity)
- Monitor usage, deprecate one if usage is minimal

---

### 4. cached_property Invalidation
**Question:** Clear cached_property on save, or remove it?

**Recommendation:** **Document as instance-scoped only**
- `cached_property` is designed for instance lifetime
- Don't fight the framework
- Add docstring: "Cached for instance lifetime only. Reload from DB to see changes."
- Consider `django-cache-memoize` if you need invalidation

---

### 5. Cloudflare Integration
**Question:** Add programmatic cache purge via Cloudflare API?

**Recommendation:** **Defer (not in scope)**
- Requires API credentials in Django settings
- Requires Cloudflare-specific code (vendor lock-in)
- HTTP caching with ETags should be sufficient
- Can add later if needed (separate project)

---

### 6. Template Fragment Caching
**Question:** Add `{% cache %}` template tags?

**Recommendation:** **Defer (not in scope)**
- Current approach (page-level HTTP caching) is simpler
- Template fragment caching adds complexity
- Django templates are fast enough for FPbase
- Consider only if profiling shows template rendering is slow

---

## Monitoring & Success Metrics

### Metrics to Track (Before/After)

**Bandwidth:**
- Average response size for `/api/proteins/spectraslugs/`
- Expected: 140 KB → <1 KB (99.7% reduction after initial load)

**Cache Hit Rates:**
- Redis cache hit rate (should increase with longer TTLs)
- Browser cache hit rate (measure via CDN/analytics)
- Expected: 50% → 80%+ hit rate

**Performance:**
- p50, p95, p99 response times for API endpoints
- Expected: Faster (fewer cache misses)

**Code Quality:**
- Lines of cache-related code (should decrease)
- Number of manual `cache.delete()` calls (should be 0)

### Monitoring Tools

**Django Debug Toolbar (local):**
- Cache panel shows cache hits/misses
- SQL panel shows query reduction

**Browser DevTools (local/production):**
- Network tab shows 304 responses
- Size column shows bandwidth savings

**Redis CLI (production):**
```bash
# Monitor cache operations
redis-cli monitor | grep -E 'GET|SET|INCR'

# Check cache hit rate
redis-cli info stats | grep keyspace_hits
```

**Sentry (production):**
- Monitor for cache-related errors
- Alert on increased cache miss rates

---

## Rollback Plan

If caching changes cause issues in production:

### Phase 2 (Spectra ETags) Rollback
```python
# Revert backend/proteins/api/views.py
@require_http_methods(["GET", "HEAD"])
def spectraslugs(request):
    spectrainfo = get_cached_spectra_info()
    response = HttpResponse(spectrainfo, content_type="application/json")
    response["Cache-Control"] = "public, max-age=600"  # No ETag
    return response
```

### Phase 3 (Signal-based invalidation) Rollback
```python
# Re-add manual cache.delete() calls
def save(self, *args, **kwargs):
    cache.delete(SPECTRA_CACHE_KEY)  # Restore manual deletion
    super().save(*args, **kwargs)
```

### Full Rollback
```bash
# Revert all caching changes
git revert <commit-range>

# Clear Redis cache (force refresh)
heroku redis:cli --app fpbase
> FLUSHDB

# Deploy
git push heroku main
```

---

## Future Improvements (Post-MVP)

1. **GraphQL Response Caching**
   - Use `django-graphql-cache` or similar
   - Per-field caching based on resolver complexity

2. **Cloudflare Cache Purge**
   - Add Cloudflare API integration
   - Purge CDN cache on version increment

3. **Stale-While-Revalidate**
   - `Cache-Control: max-age=600, stale-while-revalidate=3600`
   - Serve stale content while fetching fresh data in background

4. **Prefetch Hints**
   - Add `<link rel="prefetch">` for likely next pages
   - Preload critical API endpoints

5. **Service Worker Caching**
   - Offline-first PWA capabilities
   - Cache API responses in IndexedDB

6. **CDN Cache Analytics**
   - Track cache hit rates by endpoint
   - Identify optimization opportunities

---

## Related Issues & PRs

- **Issue #361:** Add ETag support for spectra list caching
- **Issue #360:** Enable gzip compression (✅ implemented)
- **FPBASE-5ZP:** Large HTTP payload (✅ resolved by gzip)

---

## Questions or Concerns?

**Q: Will this break existing clients?**
A: No, all changes are backward compatible. Clients that don't send `If-None-Match` get normal 200 responses.

**Q: What if Redis is flushed?**
A: Versions reset to 1, clients re-download data once, then ETags work again. No permanent impact.

**Q: How do I test ETags locally?**
A: Use `curl`:
```bash
# First request (get ETag)
curl -i http://localhost:8000/api/proteins/spectraslugs/

# Second request (conditional)
curl -i -H 'If-None-Match: "123"' http://localhost:8000/api/proteins/spectraslugs/
```

**Q: What about browser compatibility?**
A: ETags are supported by all modern browsers (IE6+, Chrome, Firefox, Safari, Edge).

---

**Last Updated:** 2025-01-09
**Status:** Planning
**Owner:** TBD
