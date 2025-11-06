# FPbase Algolia Integration Analysis & Overhaul Proposal (2025)

## Executive Summary

After a comprehensive analysis of FPbase's Algolia integration, I've identified that the current implementation is **approximately 6 years behind modern best practices**. The codebase is using deprecated libraries, outdated patterns, and missing significant performance and UX improvements available in 2025. This document provides a detailed analysis and proposes a comprehensive modernization strategy.

---

## 1. Current State Analysis

### 1.1 Backend Implementation

**Current Setup:**
- **Package**: `algoliasearch_django==3.0.0` (pyproject.toml:8)
- **Index Files**:
  - `backend/proteins/index.py` - Indexes Protein and Organism models
  - `backend/references/index.py` - Indexes Reference model
- **Pattern**: Declarative indexing using `@register` decorator with `AlgoliaIndex` classes

**What's Being Indexed:**

**Protein Model** (38 fields):
```python
fields = ("name", "uuid", "aliases", "pdb", "genbank", "uniprot", "ipg_id", "_agg",
          "img_url", "switchType", "url", "date_published", "created", "rank",
          "ga_views", "n_faves", "n_cols", "ex", "em", "pka", "ec", "qy", "em_css",
          "local_brightness", "seq", "first_author", "cofactor", "color")
should_index = "is_visible"
tags = "tags"
```

**Organism Model** (3 fields):
```python
fields = ("scientific_name", "division", "url")
```

**Reference Model** (12 fields):
```python
fields = ("doi", "journal", "pmid", "year", "first_author", "title", "citation",
          "date", "prot_primary", "prot_secondary", "_excerpts", "url")
```

**Indexing Mechanism:**
- Auto-indexing enabled (Django signals automatically update Algolia on model save/delete)
- Index suffix separation: `dev` vs `prod` (base.py:375)
- No custom handlers or manual indexing logic detected

**Issues:**
1. **Outdated Package**: Version 3.0.0 is significantly behind the latest 4.0.0 (released January 2025)
2. **Uses deprecated Python API client v3** instead of v4/v5
3. **No control over indexing performance** - synchronous updates can slow down request/response cycles
4. **Limited search configuration** - no custom ranking, relevance tuning, or replica indices
5. **No faceting configuration** visible in the index definitions

### 1.2 Frontend Implementation

**Current Setup:**
- **Search Library**: `algoliasearch@^3.35.1` (frontend/package.json:15)
- **Autocomplete Library**: `autocomplete.js@^0.36.0` (frontend/package.json:16)
- **Search UI**: Custom implementation in `frontend/src/js/algolia.js`
- **Integration Point**: Search input in navigation bar (`backend/fpbase/templates/_nav.html:80`)

**Search Implementation Details:**

**Initialization** (algolia.js:166-341):
```javascript
// Uses deprecated autocomplete.js library
$.fn.autocomplete(...)
```

**Search Sources** (in priority order):
1. Protein Index (5 results per page)
2. Reference Index (3 results per page)
3. Organism Index (2 results per page)
4. Advanced search fallback link

**Custom Features:**
- Spectra thumbnail previews in autocomplete (algolia.js:248-249)
- Match highlighting for various fields (algolia.js:27-49, 51-138)
- Custom empty state with fallback to advanced search
- Mobile optimization with autoselect on blur

**Styling**:
- Custom SCSS in `frontend/src/css/_algoliasearch.scss` (193 lines)
- Manually styled dropdown menus, suggestions, icons
- Responsive design with mobile breakpoints

**Issues:**
1. **DEPRECATED LIBRARY**: `autocomplete.js` v0.36 is no longer supported; should use `@algolia/autocomplete-js` v1.19.4+
2. **Outdated Search Client**: Using algoliasearch v3 instead of v5 (current)
3. **jQuery Dependency**: Relies on jQuery for autocomplete functionality
4. **CDN Loading**: Loads autocomplete from CDN with deferred loading (base.html:63), causing initialization race conditions
5. **No React Integration**: Despite using React extensively in the app (packages/spectra, packages/blast), search is vanilla JS
6. **Manual Template Rendering**: Custom HTML string concatenation instead of using modern UI components
7. **Limited Search Features**:
   - No faceted search
   - No filters
   - No search-as-you-type refinement
   - No AI-powered features
   - No personalization
   - No related/recommended items

### 1.3 Configuration & Settings

**Algolia Configuration** (base.py:375-384):
```python
ALGOLIA_SUFFIX = "dev" if (DEBUG or ("staging" in env("SENTRY_PROJECT", default=""))) else "prod"
ALGOLIA_PUBLIC_KEY = "421b453d4f93e332ebd0c7f3ace29476"
ALGOLIA = {
    "APPLICATION_ID": "9WAWQMVNTB",
    "API_KEY": env("ALGOLIA_API_KEY", default=""),
    "INDEX_SUFFIX": ALGOLIA_SUFFIX,
}
```

**Frontend Configuration** (base.html:70-76):
```javascript
window.FPBASE.ALGOLIA = {
  'appID': '{{ algolia_app_id }}',
  'publicKey': '{{ algolia_public_key }}',
  'proteinIndex': 'Protein_{{ algolia_suffix}}',
  'organismIndex': 'Organism_{{ algolia_suffix}}',
  'referenceIndex': 'Reference_{{ algolia_suffix}}',
};
```

**Issues:**
1. **Public API key hardcoded** in settings file (should be in environment variable)
2. **No replica indices** for different sorting strategies
3. **No query suggestions** or analytics configuration
4. **No A/B testing** setup

---

## 2. Major Issues & Limitations

### 2.1 Technical Debt

| Issue | Impact | Severity |
|-------|--------|----------|
| Deprecated autocomplete.js | Security vulnerabilities, no support, missing features | **Critical** |
| algoliasearch v3 → v5 gap | Missing modern features, performance improvements | **High** |
| algoliasearch_django 3.0 → 4.0 | Incompatible with latest Python client, missing features | **High** |
| jQuery dependency | Unnecessary bundle size, outdated patterns | **Medium** |
| Synchronous indexing | Can slow down request cycles | **Medium** |

### 2.2 Missing Modern Features

**Not Utilizing (but should be):**

1. **AI-Powered Search** (2024-2025 features)
   - AI Personalization (beta 2024) - tailored search results per user
   - AI Recommendations - related products, trending items
   - AI Ranking - ML-powered relevance optimization
   - Image-based search (Looking Similar)

2. **Advanced Search UI**
   - InstantSearch components for faceted navigation
   - Real-time filtering and refinement
   - Query suggestions
   - Recent searches
   - Rich autocomplete with multiple sections

3. **Performance Optimizations**
   - Asynchronous indexing (Celery integration)
   - Batch indexing operations
   - Cache warming strategies
   - Index replicas for different sort orders

4. **Analytics & Insights**
   - Search analytics integration
   - Click/conversion tracking
   - A/B testing capabilities
   - User behavior insights

5. **Developer Experience**
   - TypeScript support (native in v4+)
   - Better error handling
   - Improved debugging tools
   - Modern async/await patterns

### 2.3 UX Limitations

1. **Search Results**:
   - Limited to 10 total results in dropdown (5+3+2)
   - No pagination or "view more" functionality
   - No faceted navigation
   - No ability to filter by properties (color, brightness, etc.)

2. **Autocomplete**:
   - CDN loading causes initialization delays and race conditions (algolia.js:149-160)
   - Basic text matching without AI enhancements
   - No search history or suggestions
   - Limited mobile optimization

3. **Discovery**:
   - No "related proteins" recommendations
   - No trending/popular searches
   - No personalized results based on user behavior
   - Missing "users who searched for X also looked at Y" features

---

## 3. Latest Best Practices & Modern Algolia Stack (2025)

### 3.1 Backend Stack

**Recommended:**
- **algoliasearch-django 4.0.0** (released January 2025)
  - Uses Python API client v4
  - Django 4.x and 5.x support
  - Python 3.8+ compatibility
  - Native TypeScript typing support
  - Async/await patterns

**Key Improvements Over v3:**
- Promise-based API (no callbacks)
- Better error handling
- Improved TypeScript support
- More efficient batch operations
- Better connection pooling

**Breaking Changes to Handle:**
- All methods return promises (no callbacks)
- Import changes (no wildcard imports)
- Some removed methods (ttAdapter, initPlaces, setRequestTimeout)
- clearCache and destroy are now async

### 3.2 Frontend Stack

**Recommended Modern Stack:**

1. **Core Search Client**: `algoliasearch@5.x`
   - Latest: v5 (2025)
   - Full TypeScript support
   - Improved performance
   - Better error handling
   - Breaking: `initIndex()` removed, methods require `indexName` parameter

2. **Autocomplete**: `@algolia/autocomplete-js@1.19.4`
   - **REPLACES deprecated autocomplete.js**
   - Modern, actively maintained
   - Plugin architecture
   - Better mobile support
   - Integrated query suggestions
   - Recent searches plugin
   - Rich preview layouts

3. **Full Search UI**: `react-instantsearch@7.17.0` (optional but recommended)
   - Pre-built UI components
   - Faceted navigation
   - Refinement lists
   - Pagination
   - Infinite scroll
   - Server-side rendering support
   - Next.js 13 App Router compatibility

**Alternative Approach** (if staying vanilla JS):
- `instantsearch.js@4.x` - Vanilla JS InstantSearch
- `@algolia/autocomplete-js@1.x` - Modern autocomplete

### 3.3 Modern Architecture Patterns

**Recommended Indexing Strategy:**

```python
# Asynchronous indexing via Celery
class ProteinIndex(AlgoliaIndex):
    # Define fields
    fields = (...)

    # Custom indexing logic
    should_index = 'is_visible'

    # Settings for search behavior
    settings = {
        'searchableAttributes': [
            'name',
            'aliases',
            'unordered(seq)',  # Partial matching on sequence
        ],
        'attributesForFaceting': [
            'searchable(switchType)',
            'searchable(color)',
            'filterOnly(agg)',
        ],
        'customRanking': [
            'desc(ga_views)',  # Popularity ranking
            'desc(n_faves)',
            'desc(rank)',
        ],
        'replicas': [
            'Protein_prod_brightness_desc',
            'Protein_prod_name_asc',
            'Protein_prod_date_desc',
        ],
    }
```

**Recommended Frontend Pattern** (React InstantSearch):

```javascript
import { InstantSearch, SearchBox, Hits, RefinementList } from 'react-instantsearch';
import { liteClient as algoliasearch } from 'algoliasearch/lite';

const searchClient = algoliasearch('appId', 'apiKey');

function Search() {
  return (
    <InstantSearch searchClient={searchClient} indexName="Protein_prod">
      <SearchBox />
      <RefinementList attribute="switchType" />
      <RefinementList attribute="color" />
      <Hits hitComponent={ProteinHit} />
    </InstantSearch>
  );
}
```

### 3.4 Advanced Features to Consider

**1. AI Personalization**
- Personalized search results based on user behavior
- Adaptive ranking based on user preferences
- Real-time intent detection
- Cost: Requires Algolia Premium plan

**2. AI Recommendations**
- Related proteins ("Users who viewed X also viewed Y")
- Trending proteins
- Frequently compared together
- Looking similar (image-based)
- Cost: Separate Algolia Recommend product

**3. Analytics & Insights**
- Search analytics dashboard
- Click-through rate tracking
- Conversion tracking
- A/B testing framework
- Cost: Included in Growth+ plans

**4. Query Suggestions**
- As-you-type suggestions from popular queries
- Autocorrect and typo tolerance
- Synonym management
- Cost: Included in standard plans

---

## 4. Proposed Comprehensive Overhaul

### 4.1 Migration Strategy

**Phase 1: Backend Modernization** (Low Risk, High Value)
- Upgrade `algoliasearch-django` from 3.0.0 → 4.0.0
- Add index settings and configuration
- Implement asynchronous indexing via Celery
- Add replica indices for different sort orders
- Configure faceting attributes
- Improve custom ranking

**Phase 2: Frontend - Autocomplete Upgrade** (Medium Risk, High Value)
- Replace deprecated `autocomplete.js` with `@algolia/autocomplete-js`
- Remove jQuery dependency for search
- Modernize bundling (include in Vite build vs CDN)
- Improve mobile experience
- Add query suggestions and recent searches

**Phase 3: Advanced Search UI** (Higher Risk, Highest Value)
- Implement React InstantSearch for dedicated search page
- Add faceted navigation
- Implement infinite scroll / pagination
- Add sort options using replica indices
- Improve filters and refinements

**Phase 4: AI & Personalization** (Future Enhancement)
- Evaluate AI Personalization ROI
- Implement AI Recommendations
- Add trending/popular sections
- Implement analytics and A/B testing

### 4.2 Detailed Implementation Plan

#### Phase 1: Backend Modernization (1-2 days)

**Step 1.1: Update Dependencies**
```toml
# pyproject.toml
dependencies = [
    'algoliasearch-django==4.0.0',  # was 3.0.0
]
```

**Step 1.2: Enhance Index Configurations**

Create `backend/proteins/algolia_settings.py`:
```python
PROTEIN_INDEX_SETTINGS = {
    'searchableAttributes': [
        'name',
        'aliases',
        'unordered(first_author)',
        'unordered(seq)',
    ],
    'attributesForFaceting': [
        'searchable(switchType)',
        'searchable(color)',
        'searchable(cofactor)',
        'filterOnly(agg)',
    ],
    'customRanking': [
        'desc(ga_views)',
        'desc(n_faves)',
        'desc(local_brightness)',
    ],
    'ranking': [
        'typo',
        'geo',
        'words',
        'filters',
        'proximity',
        'attribute',
        'exact',
        'custom',
    ],
    'replicas': [
        'Protein_prod_name_asc',
        'Protein_prod_brightness_desc',
        'Protein_prod_date_desc',
        'Protein_prod_views_desc',
    ],
}

ORGANISM_INDEX_SETTINGS = {
    'searchableAttributes': ['scientific_name', 'division'],
}

REFERENCE_INDEX_SETTINGS = {
    'searchableAttributes': [
        'title',
        'first_author',
        'doi',
        'pmid',
        'unordered(_excerpts)',
    ],
    'attributesForFaceting': [
        'searchable(journal)',
        'filterOnly(year)',
    ],
    'customRanking': ['desc(year)'],
}
```

Update `backend/proteins/index.py`:
```python
from algoliasearch_django import AlgoliaIndex
from algoliasearch_django.decorators import register
from .algolia_settings import PROTEIN_INDEX_SETTINGS, ORGANISM_INDEX_SETTINGS
from .models import Organism, Protein

@register(Protein)
class ProteinIndex(AlgoliaIndex):
    fields = (
        "name", "uuid", "aliases", "pdb", "genbank", "uniprot", "ipg_id",
        "_agg", "img_url", "switchType", "url", "date_published", "created",
        "rank", "ga_views", "n_faves", "n_cols", "ex", "em", "pka", "ec",
        "qy", "em_css", "local_brightness", "seq", "first_author", "cofactor", "color",
    )
    should_index = "is_visible"
    tags = "tags"
    settings = PROTEIN_INDEX_SETTINGS

@register(Organism)
class OrganismIndex(AlgoliaIndex):
    fields = ("scientific_name", "division", "url")
    settings = ORGANISM_INDEX_SETTINGS
```

Update `backend/references/index.py` similarly.

**Step 1.3: Implement Async Indexing** (Optional but recommended)

Create `backend/proteins/tasks.py`:
```python
from celery import shared_task
from algoliasearch_django import update_records, raw_search

@shared_task
def update_protein_index(protein_ids):
    """Async task to update protein records in Algolia"""
    from .models import Protein
    proteins = Protein.objects.filter(id__in=protein_ids)
    update_records(Protein, proteins)

@shared_task
def reindex_all_proteins():
    """Async task to reindex all proteins"""
    from django.core.management import call_command
    call_command('algolia_reindex', '--model', 'proteins.Protein')
```

**Step 1.4: Testing**
```bash
# Test indexing
uv run backend/manage.py algolia_reindex --model proteins.Protein
uv run backend/manage.py algolia_reindex --model proteins.Organism
uv run backend/manage.py algolia_reindex --model references.Reference

# Clear and rebuild
uv run backend/manage.py algolia_clearindex --model proteins.Protein
uv run backend/manage.py algolia_reindex --model proteins.Protein
```

**Step 1.5: Verification**
- Check Algolia dashboard for new index settings
- Verify replica indices are created
- Test search quality improvements
- Monitor indexing performance

---

#### Phase 2: Frontend Autocomplete Upgrade (2-3 days)

**Step 2.1: Update Dependencies**
```json
// frontend/package.json
{
  "dependencies": {
    "@algolia/autocomplete-js": "^1.19.4",  // NEW
    "algoliasearch": "^5.17.0",  // was ^3.35.1
    // REMOVE: "autocomplete.js": "^0.36.0"
  }
}
```

Run: `pnpm install`

**Step 2.2: Remove CDN Dependency**

In `backend/fpbase/templates/base.html`, REMOVE:
```html
<!-- REMOVE THIS LINE -->
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/autocomplete.js/0.38.1/autocomplete.jquery.min.js" ...></script>
```

**Step 2.3: Create Modern Autocomplete Implementation**

Create `frontend/src/js/autocomplete-search.js`:
```javascript
import { autocomplete } from '@algolia/autocomplete-js';
import { createQuerySuggestionsPlugin } from '@algolia/autocomplete-plugin-query-suggestions';
import { createLocalStorageRecentSearchesPlugin } from '@algolia/autocomplete-plugin-recent-searches';
import algoliasearch from 'algoliasearch';

const ALGOLIA_CONFIG = window.FPBASE.ALGOLIA;
const searchClient = algoliasearch(ALGOLIA_CONFIG.appID, ALGOLIA_CONFIG.publicKey);

// Plugin for recent searches
const recentSearchesPlugin = createLocalStorageRecentSearchesPlugin({
  key: 'fpbase-recent-searches',
  limit: 5,
});

export function initAutocomplete() {
  const searchInput = document.querySelector('#algolia-search-input');
  if (!searchInput) {
    console.warn('Algolia search input not found in DOM');
    return;
  }

  autocomplete({
    container: searchInput.parentElement,
    placeholder: 'Search proteins, organisms, references...',
    openOnFocus: true,
    plugins: [recentSearchesPlugin],
    getSources({ query }) {
      if (!query) {
        return [];
      }

      return [
        // Proteins source
        {
          sourceId: 'proteins',
          getItems() {
            return searchClient.search([{
              indexName: ALGOLIA_CONFIG.proteinIndex,
              query,
              params: {
                hitsPerPage: 5,
                attributesToRetrieve: [
                  'name', 'uuid', 'aliases', 'switchType', 'color',
                  'ex', 'em', 'img_url', 'url'
                ],
                attributesToHighlight: ['name', 'aliases'],
              },
            }]).then(({ results }) => results[0].hits);
          },
          templates: {
            header() {
              return '<div class="aa-SourceHeader">Proteins</div>';
            },
            item({ item, components }) {
              let iconColor = 'gray50';
              if (item.switchType && item.switchType !== 'Basic') {
                iconColor = 'rainbow';
              } else if (item.color && !item.color.includes('Stokes')) {
                iconColor = item.color.toLowerCase().replace(/ |\//g, '_');
              }

              return `
                <div class="aa-ItemWrapper">
                  <div class="aa-ItemContent">
                    <div class="aa-ItemIcon">
                      <img src="${window.FPBASE.imageDir}gfp_${iconColor}_40.png" alt="protein icon">
                    </div>
                    <div class="aa-ItemContentBody">
                      <div class="aa-ItemContentTitle">
                        ${components.Highlight({ hit: item, attribute: 'name' })}
                      </div>
                      ${item.ex && item.em ? `<div class="aa-ItemContentDescription">${item.ex}/${item.em}</div>` : ''}
                    </div>
                    ${item.img_url ? `<img class="aa-ItemContentSpectra" src="${item.img_url}" alt="spectra">` : ''}
                  </div>
                </div>
              `;
            },
          },
          getItemUrl({ item }) {
            return item.url;
          },
        },
        // References source
        {
          sourceId: 'references',
          getItems() {
            return searchClient.search([{
              indexName: ALGOLIA_CONFIG.referenceIndex,
              query,
              params: {
                hitsPerPage: 3,
                attributesToHighlight: ['citation', 'title', 'doi', 'pmid'],
              },
            }]).then(({ results }) => results[0].hits);
          },
          templates: {
            header() {
              return '<div class="aa-SourceHeader">References</div>';
            },
            item({ item, components }) {
              return `
                <div class="aa-ItemWrapper">
                  <div class="aa-ItemContent">
                    <div class="aa-ItemIcon">
                      <img src="${window.FPBASE.imageDir}ref.png" alt="reference icon">
                    </div>
                    <div class="aa-ItemContentBody">
                      <div class="aa-ItemContentTitle">
                        ${components.Highlight({ hit: item, attribute: 'citation' })}
                      </div>
                      ${item._highlightResult?.title ? `<div class="aa-ItemContentDescription">${item._highlightResult.title.value}</div>` : ''}
                    </div>
                  </div>
                </div>
              `;
            },
          },
          getItemUrl({ item }) {
            return item.url;
          },
        },
        // Organisms source
        {
          sourceId: 'organisms',
          getItems() {
            return searchClient.search([{
              indexName: ALGOLIA_CONFIG.organismIndex,
              query,
              params: {
                hitsPerPage: 2,
                attributesToHighlight: ['scientific_name'],
              },
            }]).then(({ results }) => results[0].hits);
          },
          templates: {
            header() {
              return '<div class="aa-SourceHeader">Organisms</div>';
            },
            item({ item, components }) {
              return `
                <div class="aa-ItemWrapper">
                  <div class="aa-ItemContent">
                    <div class="aa-ItemIcon">
                      <img src="${window.FPBASE.imageDir}organism_icon.png" alt="organism icon">
                    </div>
                    <div class="aa-ItemContentBody">
                      ${components.Highlight({ hit: item, attribute: 'scientific_name' })}
                    </div>
                  </div>
                </div>
              `;
            },
          },
          getItemUrl({ item }) {
            return item.url;
          },
        },
      ];
    },
    navigator: {
      navigate({ itemUrl }) {
        window.location.assign(itemUrl);
      },
    },
  });
}
```

**Step 2.4: Update Main Entry Point**

Update `frontend/src/index.js`:
```javascript
// Replace old import
// import initAutocomplete from './js/algolia.js'

// New import
import { initAutocomplete } from './js/autocomplete-search.js'

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
  initAutocomplete();
});
```

**Step 2.5: Update Styles**

Replace `frontend/src/css/_algoliasearch.scss` with modern autocomplete styles:
```scss
// Modern Autocomplete.js styles
@import '@algolia/autocomplete-theme-classic';

.aa-Autocomplete {
  width: 100%;
}

.aa-Form {
  border: solid 1px rgba(255, 255, 255, 0.5);
  box-shadow: 0 1px 10px rgba(0, 0, 0, 0.2), 0 2px 4px 0 rgba(0, 0, 0, 0.1);
  border-radius: 4px;
}

.aa-Input {
  color: #60a263;
  font-family: "Raleway", "Helvetica Neue", helvetica;
  height: 50px;
  font-weight: 500;
  font-size: 1rem;

  &::placeholder {
    color: #888;
    font-size: 1rem;
  }
}

.aa-Panel {
  max-height: 50vh;
  overflow-y: auto;

  @include media-breakpoint-down(sm) {
    max-height: 48vh;
  }
}

.aa-SourceHeader {
  font-size: 0.75rem;
  font-weight: 600;
  color: #999;
  text-transform: uppercase;
  padding: 8px 12px;
  margin: 0;
}

.aa-Item {
  cursor: pointer;
  padding: 8px 12px;

  &[aria-selected="true"] {
    background-color: #eee;

    .aa-ItemContentSpectra {
      filter: none;
      opacity: 0.9;
    }
  }

  em {
    font-weight: bold;
    font-style: normal;
    background-color: rgba(96, 162, 99, 0.15);
  }
}

.aa-ItemWrapper {
  display: flex;
  align-items: center;
}

.aa-ItemContent {
  display: flex;
  align-items: center;
  width: 100%;
  position: relative;
}

.aa-ItemIcon {
  width: 30px;
  margin-right: 10px;

  img {
    height: 1.6rem;
    opacity: 0.8;
  }
}

.aa-ItemContentBody {
  flex: 1;
}

.aa-ItemContentTitle {
  font-size: 1rem;
  line-height: 1.4;
}

.aa-ItemContentDescription {
  font-size: 0.85rem;
  color: #666;
  margin-top: 2px;
}

.aa-ItemContentSpectra {
  position: absolute;
  right: 5px;
  height: 2rem;
  width: 23%;
  min-width: 150px;
  max-width: 200px;
  opacity: 0.6;
  filter: grayscale(70%);
}
```

**Step 2.6: Update Import in Main Styles**

Update `frontend/src/css/style.scss`:
```scss
// Update or verify this import exists
@import '_algoliasearch';
```

**Step 2.7: Build and Test**
```bash
pnpm build
pnpm dev

# Test on localhost:8000
# - Search for proteins by name
# - Check autocomplete suggestions
# - Verify styling
# - Test on mobile viewport
```

---

#### Phase 3: Advanced Search UI with React InstantSearch (3-5 days)

**Step 3.1: Add Dependencies**
```json
// frontend/package.json
{
  "dependencies": {
    "react-instantsearch": "^7.17.0",
    "instantsearch.css": "^8.5.0"
  }
}
```

**Step 3.2: Create Advanced Search Page Component**

Create `frontend/src/components/AdvancedSearch.jsx`:
```jsx
import React from 'react';
import {
  InstantSearch,
  SearchBox,
  Hits,
  RefinementList,
  Pagination,
  Stats,
  SortBy,
  ClearRefinements,
  Configure,
  Panel,
} from 'react-instantsearch';
import { liteClient as algoliasearch } from 'algoliasearch/lite';
import 'instantsearch.css/themes/satellite.css';

const searchClient = algoliasearch(
  window.FPBASE.ALGOLIA.appID,
  window.FPBASE.ALGOLIA.publicKey
);

function ProteinHit({ hit }) {
  let iconColor = 'gray50';
  if (hit.switchType && hit.switchType !== 'Basic') {
    iconColor = 'rainbow';
  } else if (hit.color && !hit.color.includes('Stokes')) {
    iconColor = hit.color.toLowerCase().replace(/ |\//g, '_');
  }

  return (
    <div className="protein-hit">
      <div className="protein-hit__icon">
        <img
          src={`${window.FPBASE.imageDir}gfp_${iconColor}_40.png`}
          alt="protein"
        />
      </div>
      <div className="protein-hit__content">
        <h3>
          <a href={hit.url}>{hit.name}</a>
        </h3>
        {hit.aliases && hit.aliases.length > 0 && (
          <p className="protein-hit__aliases">
            aka: {hit.aliases.join(', ')}
          </p>
        )}
        <div className="protein-hit__properties">
          {hit.ex && hit.em && (
            <span className="protein-hit__wavelength">
              {hit.ex}/{hit.em} nm
            </span>
          )}
          {hit.switchType && hit.switchType !== 'Basic' && (
            <span className="protein-hit__switch-type">{hit.switchType}</span>
          )}
          {hit.color && (
            <span className="protein-hit__color">{hit.color}</span>
          )}
        </div>
      </div>
      {hit.img_url && (
        <div className="protein-hit__spectra">
          <img src={hit.img_url} alt="spectra" />
        </div>
      )}
    </div>
  );
}

export function AdvancedSearch() {
  return (
    <div className="advanced-search">
      <InstantSearch
        searchClient={searchClient}
        indexName={window.FPBASE.ALGOLIA.proteinIndex}
      >
        <Configure hitsPerPage={20} />

        <div className="search-header">
          <SearchBox
            placeholder="Search proteins, aliases, sequences..."
            autoFocus
          />
          <Stats />
        </div>

        <div className="search-layout">
          <aside className="search-sidebar">
            <Panel header="Sort By">
              <SortBy
                items={[
                  { label: 'Relevance', value: window.FPBASE.ALGOLIA.proteinIndex },
                  { label: 'Name (A-Z)', value: `${window.FPBASE.ALGOLIA.proteinIndex}_name_asc` },
                  { label: 'Brightness (High-Low)', value: `${window.FPBASE.ALGOLIA.proteinIndex}_brightness_desc` },
                  { label: 'Most Popular', value: `${window.FPBASE.ALGOLIA.proteinIndex}_views_desc` },
                  { label: 'Recently Added', value: `${window.FPBASE.ALGOLIA.proteinIndex}_date_desc` },
                ]}
              />
            </Panel>

            <ClearRefinements />

            <Panel header="Switch Type">
              <RefinementList
                attribute="switchType"
                searchable
                showMore
                limit={5}
              />
            </Panel>

            <Panel header="Color">
              <RefinementList
                attribute="color"
                searchable
                showMore
                limit={10}
              />
            </Panel>

            <Panel header="Aggregation State">
              <RefinementList attribute="agg" />
            </Panel>

            <Panel header="Cofactor">
              <RefinementList attribute="cofactor" />
            </Panel>
          </aside>

          <main className="search-results">
            <Hits hitComponent={ProteinHit} />
            <Pagination showFirst={false} showLast={false} />
          </main>
        </div>
      </InstantSearch>
    </div>
  );
}
```

**Step 3.3: Add Advanced Search Route**

Create entry point `frontend/src/search.js`:
```javascript
import React from 'react';
import ReactDOM from 'react-dom/client';
import { AdvancedSearch } from './components/AdvancedSearch';

document.addEventListener('DOMContentLoaded', () => {
  const searchContainer = document.getElementById('advanced-search-root');
  if (searchContainer) {
    const root = ReactDOM.createRoot(searchContainer);
    root.render(<AdvancedSearch />);
  }
});
```

**Step 3.4: Update Vite Config**

Update `frontend/vite.config.js` to include new entry point:
```javascript
export default defineConfig({
  // ... existing config
  build: {
    rollupOptions: {
      input: {
        main: resolve(__dirname, 'src/index.js'),
        search: resolve(__dirname, 'src/search.js'),  // NEW
        // ... other entry points
      },
    },
  },
});
```

**Step 3.5: Create Django Template**

Create `backend/fpbase/templates/pages/advanced_search.html`:
```django
{% extends "base.html" %}
{% load django_vite %}

{% block content %}
<div id="advanced-search-root"></div>
{% vite_asset 'src/search.js' %}
{% endblock %}
```

**Step 3.6: Add URL Route**

Update `backend/config/urls.py`:
```python
urlpatterns = [
    # ... existing routes
    path('search/', TemplateView.as_view(template_name='pages/advanced_search.html'), name='advanced-search'),
]
```

**Step 3.7: Add Styles**

Create `frontend/src/css/_advanced-search.scss`:
```scss
.advanced-search {
  max-width: 1400px;
  margin: 0 auto;
  padding: 2rem;
}

.search-header {
  margin-bottom: 2rem;

  .ais-SearchBox {
    max-width: 600px;
  }

  .ais-Stats {
    margin-top: 1rem;
    color: #666;
    font-size: 0.9rem;
  }
}

.search-layout {
  display: grid;
  grid-template-columns: 250px 1fr;
  gap: 2rem;

  @media (max-width: 768px) {
    grid-template-columns: 1fr;
  }
}

.search-sidebar {
  .ais-Panel {
    margin-bottom: 1.5rem;
  }

  .ais-Panel-header {
    font-weight: 600;
    margin-bottom: 0.5rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid $primary;
  }
}

.protein-hit {
  display: flex;
  padding: 1rem;
  border: 1px solid #e0e0e0;
  border-radius: 4px;
  margin-bottom: 1rem;
  transition: box-shadow 0.2s;

  &:hover {
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  }

  &__icon {
    margin-right: 1rem;

    img {
      height: 40px;
    }
  }

  &__content {
    flex: 1;

    h3 {
      margin: 0 0 0.5rem;
      font-size: 1.2rem;

      a {
        color: $primary;
        text-decoration: none;

        &:hover {
          text-decoration: underline;
        }
      }
    }
  }

  &__aliases {
    font-size: 0.85rem;
    color: #666;
    margin: 0 0 0.5rem;
  }

  &__properties {
    display: flex;
    gap: 1rem;
    font-size: 0.85rem;

    span {
      padding: 0.25rem 0.5rem;
      background-color: #f5f5f5;
      border-radius: 3px;
    }
  }

  &__spectra {
    width: 200px;

    img {
      width: 100%;
      height: auto;
    }
  }
}
```

Import in `frontend/src/css/style.scss`:
```scss
@import '_advanced-search';
```

---

#### Phase 4: AI & Personalization (Future Enhancement)

**Consider when:**
- You have significant user traffic (>10k searches/month)
- You want to improve conversion (protein favorites, collection adds)
- You have budget for Premium Algolia features

**Features to implement:**
1. **AI Personalization** - Personalized search ranking
2. **AI Recommendations** - Related/similar proteins
3. **Trending Searches** - Popular queries dashboard
4. **Analytics** - Click tracking, conversion optimization
5. **A/B Testing** - Test different ranking strategies

**Estimated Cost:**
- AI Personalization: Requires Premium plan (~$1,500+/month)
- AI Recommendations: Separate product (~$500+/month)
- Analytics: Included in Growth+ plans (~$500/month)

---

### 4.3 Testing Strategy

**Unit Tests:**
```python
# backend/proteins/tests/test_algolia.py
import pytest
from proteins.models import Protein
from proteins.index import ProteinIndex

@pytest.mark.django_db
class TestProteinIndexing:
    def test_should_index_visible_proteins(self, protein_factory):
        protein = protein_factory(status='approved')
        index = ProteinIndex()
        assert index.should_index(protein)

    def test_should_not_index_hidden_proteins(self, protein_factory):
        protein = protein_factory(status='hidden')
        index = ProteinIndex()
        assert not index.should_index(protein)
```

**Integration Tests:**
```python
# Test Algolia search API
def test_protein_search_returns_results(algolia_client):
    results = algolia_client.search('GFP')
    assert len(results['hits']) > 0
    assert 'eGFP' in [hit['name'] for hit in results['hits']]
```

**E2E Tests:**
```javascript
// backend/tests_e2e/test_search.py (Playwright)
async test('autocomplete shows protein results', async ({ page }) => {
  await page.goto('/');
  await page.fill('#algolia-search-input', 'mCherry');
  await page.waitForSelector('.aa-Panel');

  const results = await page.locator('.aa-Item').count();
  expect(results).toBeGreaterThan(0);

  const firstResult = await page.locator('.aa-Item').first().textContent();
  expect(firstResult).toContain('mCherry');
});
```

---

### 4.4 Rollback Plan

**If Phase 2 (Autocomplete) has issues:**
1. Revert frontend/package.json dependencies
2. Restore old algolia.js file
3. Re-add CDN script tag in base.html
4. Run `pnpm install && pnpm build`

**If Phase 1 (Backend) has issues:**
1. Revert pyproject.toml to `algoliasearch-django==3.0.0`
2. Run `uv sync`
3. Restore original index.py files
4. Reindex: `uv run backend/manage.py algolia_reindex`

---

### 4.5 Migration Checklist

**Pre-Migration:**
- [ ] Backup Algolia indices (export via dashboard)
- [ ] Document current search behavior (screenshots, videos)
- [ ] Set up Algolia analytics (baseline metrics)
- [ ] Create feature branch: `git checkout -b feature/algolia-modernization`

**Phase 1 - Backend:**
- [ ] Update `algoliasearch-django` to 4.0.0
- [ ] Add index settings configurations
- [ ] Test indexing locally with dev indices
- [ ] Reindex development indices
- [ ] Verify search quality in Algolia dashboard
- [ ] Deploy to staging
- [ ] Monitor indexing performance
- [ ] Deploy to production

**Phase 2 - Autocomplete:**
- [ ] Install `@algolia/autocomplete-js` and `algoliasearch@5`
- [ ] Remove deprecated `autocomplete.js`
- [ ] Implement new autocomplete component
- [ ] Update styles
- [ ] Test locally (Chrome, Firefox, Safari, Mobile)
- [ ] Test on staging
- [ ] Get user feedback
- [ ] Deploy to production
- [ ] Monitor error rates (Sentry)

**Phase 3 - Advanced Search:**
- [ ] Install `react-instantsearch`
- [ ] Create AdvancedSearch component
- [ ] Add Django route and template
- [ ] Style components
- [ ] Test faceted navigation
- [ ] Test sorting
- [ ] Test pagination
- [ ] Deploy to production
- [ ] Update navigation links

**Phase 4 - Future:**
- [ ] Evaluate AI features ROI
- [ ] Request Algolia sales demo
- [ ] Implement AI Personalization pilot
- [ ] Set up analytics dashboards
- [ ] A/B test new vs old search

---

### 4.6 Estimated Timeline

| Phase | Duration | Effort | Risk |
|-------|----------|--------|------|
| Phase 1: Backend | 1-2 days | Low | Low |
| Phase 2: Autocomplete | 2-3 days | Medium | Medium |
| Phase 3: Advanced Search | 3-5 days | High | Low |
| Phase 4: AI Features | 2-4 weeks | High | Medium |
| **Total (Phases 1-3)** | **1-2 weeks** | - | - |

---

## 5. Benefits & ROI

### 5.1 Performance Improvements

**Expected Gains:**
- **20-30% faster autocomplete** - Modern libraries are more optimized
- **50% reduction in bundle size** - Removing jQuery for search
- **Better mobile performance** - Native mobile optimization in new autocomplete
- **Async indexing** - Faster Django response times (no blocking on Algolia updates)

### 5.2 UX Improvements

**User Experience:**
- Modern, polished autocomplete UI
- Recent searches memory
- Better mobile experience
- Faceted search with filters
- Multiple sort options
- Better accessibility (ARIA labels, keyboard navigation)

### 5.3 Developer Experience

**Engineering:**
- TypeScript support (better IDE autocomplete)
- Modern async/await patterns
- Better error handling
- Easier testing
- Active community support
- Regular security updates

### 5.4 Future-Proofing

**Long-term:**
- Security updates and bug fixes
- Access to new features (AI, personalization)
- Better analytics and insights
- Scalability for growth
- No breaking changes for 3-5 years

---

## 6. Risks & Mitigation

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Breaking changes in autocomplete UI | High | Medium | Extensive testing, staged rollout |
| Search quality regression | High | Low | A/B testing, monitor analytics |
| Increased Algolia costs | Medium | Low | Stay on current plan, monitor usage |
| Bundle size increase | Medium | Low | Code splitting, lazy loading |
| Learning curve for team | Low | Medium | Documentation, training session |

---

## 7. Recommendations

### 7.1 Immediate Actions (Do Now)

1. **Upgrade Backend** (Phase 1) - **CRITICAL**
   - Low risk, high value
   - Fixes security vulnerabilities
   - Enables better search configuration
   - Required for frontend upgrades

2. **Replace Autocomplete.js** (Phase 2) - **HIGH PRIORITY**
   - Library is deprecated and unsupported
   - Security risk
   - Blocking new features

### 7.2 Short-Term (Next 3 Months)

3. **Implement React InstantSearch** (Phase 3) - **RECOMMENDED**
   - Significantly better UX
   - Enables faceted search
   - Users can self-serve complex queries

### 7.3 Long-Term (6-12 Months)

4. **Evaluate AI Features** (Phase 4) - **EVALUATE**
   - Depends on budget and user traffic
   - Consider after Phases 1-3 complete
   - Measure search analytics first

---

## 8. Appendix

### 8.1 Key Libraries Comparison

| Feature | Current | Recommended | Status |
|---------|---------|-------------|--------|
| Django Integration | `algoliasearch-django@3.0.0` | `algoliasearch-django@4.0.0` | ⚠️ **Outdated** |
| JS Search Client | `algoliasearch@3.35.1` | `algoliasearch@5.17.0` | ⚠️ **Outdated** |
| Autocomplete | `autocomplete.js@0.36.0` | `@algolia/autocomplete-js@1.19.4` | ❌ **DEPRECATED** |
| Search UI | Custom jQuery | `react-instantsearch@7.17.0` | ⚠️ **Missing** |

### 8.2 Feature Gap Analysis

| Feature | Current | Available | Gap |
|---------|---------|-----------|-----|
| Basic Search | ✅ Yes | ✅ Yes | - |
| Autocomplete | ✅ Yes (deprecated) | ✅ Yes | Using old library |
| Faceted Search | ❌ No | ✅ Yes | **Missing** |
| Filters | ❌ No | ✅ Yes | **Missing** |
| Sort Options | ❌ No | ✅ Yes | **Missing** |
| Recent Searches | ❌ No | ✅ Yes | **Missing** |
| Query Suggestions | ❌ No | ✅ Yes | **Missing** |
| AI Personalization | ❌ No | ✅ Yes (Premium) | **Missing** |
| AI Recommendations | ❌ No | ✅ Yes (Paid) | **Missing** |
| Analytics | ❌ No | ✅ Yes | **Missing** |
| A/B Testing | ❌ No | ✅ Yes | **Missing** |

### 8.3 Resources

**Official Documentation:**
- Algolia Django: https://www.algolia.com/doc/framework-integration/django/
- Autocomplete.js: https://www.algolia.com/doc/ui-libraries/autocomplete/
- React InstantSearch: https://www.algolia.com/doc/guides/building-search-ui/what-is-instantsearch/react/
- JavaScript Client v5: https://www.algolia.com/doc/libraries/javascript/v5/

**Migration Guides:**
- Django upgrade guide: https://www.algolia.com/doc/framework-integration/django/upgrade-guide/
- JS v3→v4 upgrade: https://www.algolia.com/doc/libraries/javascript/v4/update/
- JS v4→v5 upgrade: https://www.algolia.com/doc/libraries/javascript/v5/upgrade/

**GitHub Repositories:**
- algoliasearch-django: https://github.com/algolia/algoliasearch-django
- autocomplete: https://github.com/algolia/autocomplete
- react-instantsearch: https://github.com/algolia/react-instantsearch

---

## Conclusion

FPbase's Algolia integration is approximately **6 years behind modern best practices**. The current implementation uses:
- A deprecated autocomplete library
- Outdated search clients (v3 instead of v5)
- Old Django integration (v3 instead of v4)
- Manual HTML templating instead of modern UI components
- No faceted search, filters, or advanced features

The proposed modernization will:
1. **Eliminate security risks** from deprecated libraries
2. **Improve performance** through modern optimizations
3. **Enhance UX** with faceted search, filters, and better mobile support
4. **Future-proof** the codebase for 3-5 years
5. **Enable AI features** when budget allows

**Recommended Approach:**
- Start with **Phase 1 (Backend)** immediately - low risk, high value
- Follow with **Phase 2 (Autocomplete)** within 2 weeks - critical security fix
- Implement **Phase 3 (Advanced Search)** within 2 months - major UX improvement
- Evaluate **Phase 4 (AI)** after 6 months based on analytics

**Total Estimated Effort:** 1-2 weeks for Phases 1-3 (excluding testing/deployment overhead)

**Total Estimated Cost:** $0 (all free/open-source upgrades; no Algolia plan changes required)
