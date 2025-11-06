# FPbase Algolia Integration: Complete From-Scratch Redesign

> **TL;DR**: Remove the stale `algoliasearch-django` wrapper entirely. Build a modern, frontend-first architecture with minimal backend indexing code using the direct Python client, React InstantSearch UI, and proper separation of concerns.

---

## Why Redesign From Scratch?

### The Problem with Current Approach

1. **algoliasearch-django is INACTIVE**
   - Last release: 7 months ago
   - Classification: "Discontinued project"
   - Only 1 dependent project
   - 19 open issues
   - Not keeping up with Python client updates

2. **Wrong Architecture Pattern**
   - Django wrapper forces backend-centric thinking
   - Auto-indexing via signals couples search to Django ORM
   - Synchronous operations block request/response cycles
   - Over-abstraction for minimal value

3. **Fighting the Framework**
   - You have React in your stack but use jQuery for search
   - Mixing template-based and SPA patterns
   - CDN-loaded deprecated libraries

### Algolia's Official Guidance (2025)

> **"Algolia highly recommends front-end search"**
> - Up to **10x faster** for end-users
> - Distributed API with **global datacenters**
> - Better **availability** (no single point of failure)
> - Only use backend when you need pre-processing, security filters, or SEO

---

## The Modern Architecture: Frontend-First

### Core Philosophy

1. **Frontend handles ALL search** → Direct to Algolia from browser
2. **Backend handles ONLY indexing** → Simple, explicit, async
3. **React for UI** → Leveraging InstantSearch components
4. **Python client directly** → No Django wrapper abstraction
5. **Celery for async** → Index updates don't block requests

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                         BROWSER                              │
│  ┌────────────────────────────────────────────────────┐     │
│  │  React InstantSearch UI                            │     │
│  │  - Autocomplete                                    │     │
│  │  - Search page with facets                        │     │
│  │  - Filters, sorting, pagination                   │     │
│  └────────────────────────────────────────────────────┘     │
│                           │                                  │
│                           │ (Direct API calls)               │
│                           ↓                                  │
│                    ┌─────────────┐                          │
│                    │   Algolia   │                          │
│                    │  (CDN-like) │                          │
│                    └─────────────┘                          │
└─────────────────────────────────────────────────────────────┘
                               ↑
                               │ (Index updates only)
                               │
┌──────────────────────────────┴──────────────────────────────┐
│                      DJANGO BACKEND                          │
│  ┌──────────────────────────────────────────────────────┐  │
│  │  Custom Indexing Service                             │  │
│  │  - proteins/algolia.py                               │  │
│  │  - Explicit index() methods                          │  │
│  │  - Batch operations                                  │  │
│  │  - Direct algoliasearch client                       │  │
│  └──────────────────────────────────────────────────────┘  │
│                           │                                  │
│                           ↓                                  │
│  ┌──────────────────────────────────────────────────────┐  │
│  │  Celery Tasks                                        │  │
│  │  - Async indexing                                    │  │
│  │  - Batch reindexing                                  │  │
│  │  - No blocking on save()                             │  │
│  └──────────────────────────────────────────────────────┘  │
│                                                              │
│  Triggered by:                                               │
│  - Model save (via custom signal handler)                   │
│  - Management commands (reindex)                             │
│  - Admin actions (bulk reindex)                              │
└──────────────────────────────────────────────────────────────┘
```

---

## Implementation Plan

### Phase 1: Backend Indexing Service (Replace algoliasearch-django)

**Goal**: Simple, explicit indexing without the Django wrapper

#### Step 1.1: Dependencies

**REMOVE:**
```toml
# pyproject.toml - DELETE THIS
'algoliasearch_django==3.0.0',
```

**ADD:**
```toml
# pyproject.toml
dependencies = [
    'algoliasearch>=4.0,<5.0',  # Direct client, actively maintained
]
```

#### Step 1.2: Create Custom Indexing Service

Create `backend/proteins/algolia.py`:

```python
"""
Custom Algolia indexing service for FPbase.

Why not algoliasearch-django?
- That package is inactive/discontinued
- Frontend-first architecture means minimal backend indexing needs
- Direct client gives us more control and is actively maintained
"""

from typing import Any, Iterable
from django.conf import settings
from algoliasearch.search.client import SearchClient


class AlgoliaIndexer:
    """Base indexing service for Algolia."""

    def __init__(self):
        self.client = SearchClient.create(
            settings.ALGOLIA['APPLICATION_ID'],
            settings.ALGOLIA['API_KEY']
        )
        self.index_suffix = settings.ALGOLIA['INDEX_SUFFIX']

    def get_index_name(self, base_name: str) -> str:
        """Get full index name with environment suffix."""
        return f"{base_name}_{self.index_suffix}"

    def get_index(self, base_name: str):
        """Get index instance."""
        return self.client.init_index(self.get_index_name(base_name))


class ProteinIndexer(AlgoliaIndexer):
    """Algolia indexing for Protein model."""

    INDEX_NAME = "Protein"

    def serialize_protein(self, protein) -> dict[str, Any]:
        """Convert protein to Algolia record."""
        from proteins.models import Protein

        if not protein.is_visible():
            return None  # Don't index hidden proteins

        # Build the record
        record = {
            'objectID': str(protein.uuid),  # Algolia requires objectID
            'name': protein.name,
            'uuid': str(protein.uuid),
            'url': protein.get_absolute_url(),
            'slug': protein.slug,

            # Searchable text
            'aliases': list(protein.aliases) if protein.aliases else [],
            'seq': str(protein.seq) if protein.seq else None,
            'first_author': protein.primary_reference.first_author if protein.primary_reference else None,

            # External IDs
            'pdb': protein.pdb,
            'genbank': protein.genbank,
            'uniprot': protein.uniprot,
            'ipg_id': protein.ipg_id,

            # Properties for faceting
            'switchType': protein.get_switch_type_display(),
            'color': protein.get_color_display() if protein.color else None,
            'cofactor': protein.cofactor,
            'agg': protein.agg,

            # Numeric attributes
            'ex': protein.default_state.ex_max if protein.default_state else None,
            'em': protein.default_state.em_max if protein.default_state else None,
            'pka': protein.pka,
            'ec': protein.default_state.ext_coeff if protein.default_state else None,
            'qy': protein.default_state.qy if protein.default_state else None,
            'local_brightness': protein.local_brightness,

            # Popularity metrics for custom ranking
            'ga_views': protein.ga_views or 0,
            'n_faves': protein.n_faves or 0,
            'n_cols': protein.n_cols or 0,
            'rank': protein.rank or 0,

            # Metadata
            'created': int(protein.created.timestamp()) if protein.created else None,
            'date_published': int(protein.date_published.timestamp()) if protein.date_published else None,

            # UI helpers
            'img_url': protein.spectra_img_url if hasattr(protein, 'spectra_img_url') else None,
            'em_css': protein.em_color_css if hasattr(protein, 'em_color_css') else None,

            # Tags for filtering
            '_tags': self._get_tags(protein),
        }

        return {k: v for k, v in record.items() if v is not None}

    def _get_tags(self, protein) -> list[str]:
        """Get tags for filtering."""
        tags = []
        if protein.switchType != 'b':  # Not basic
            tags.append('switchable')
        if protein.agg == 'm':
            tags.append('monomer')
        if protein.color:
            tags.append(f'color:{protein.color}')
        return tags

    def index_protein(self, protein) -> None:
        """Index a single protein."""
        record = self.serialize_protein(protein)
        if record:
            index = self.get_index(self.INDEX_NAME)
            index.save_object(record)
        else:
            # Protein shouldn't be indexed (hidden), delete if exists
            self.delete_protein(protein)

    def index_proteins_batch(self, proteins: Iterable) -> None:
        """Index multiple proteins in batch."""
        records = []
        for protein in proteins:
            record = self.serialize_protein(protein)
            if record:
                records.append(record)

        if records:
            index = self.get_index(self.INDEX_NAME)
            index.save_objects(records)

    def delete_protein(self, protein) -> None:
        """Remove protein from index."""
        index = self.get_index(self.INDEX_NAME)
        index.delete_object(str(protein.uuid))

    def clear_index(self) -> None:
        """Clear all records from index."""
        index = self.get_index(self.INDEX_NAME)
        index.clear_objects()

    def configure_index(self) -> None:
        """Configure index settings (searchable attributes, ranking, etc)."""
        index = self.get_index(self.INDEX_NAME)

        index.set_settings({
            # What fields to search in (order matters for relevance)
            'searchableAttributes': [
                'name',           # Highest priority
                'aliases',
                'unordered(first_author)',
                'unordered(seq)',  # Unordered = position doesn't affect ranking
            ],

            # What can be used for faceting/filtering
            'attributesForFaceting': [
                'searchable(switchType)',  # Also searchable within facets
                'searchable(color)',
                'searchable(cofactor)',
                'filterOnly(agg)',  # Can filter but not search
                'filterOnly(_tags)',
            ],

            # Custom ranking criteria (applied after textual relevance)
            'customRanking': [
                'desc(ga_views)',        # Most viewed first
                'desc(n_faves)',         # Most favorited
                'desc(local_brightness)', # Brightest
                'desc(rank)',
            ],

            # Faceting display order
            'renderingContent': {
                'facetOrdering': {
                    'facets': {
                        'order': ['switchType', 'color', 'cofactor', 'agg']
                    }
                }
            },

            # Enable de-duplication by name
            'attributeForDistinct': 'name',
            'distinct': 1,

            # Pagination
            'hitsPerPage': 20,
            'paginationLimitedTo': 1000,

            # Typo tolerance
            'typoTolerance': True,
            'minWordSizefor1Typo': 4,
            'minWordSizefor2Typos': 8,
        })

        # Create replica indices for different sort orders
        self._create_replicas()

    def _create_replicas(self) -> None:
        """Create replica indices for alternative sorting."""
        base_index = self.get_index(self.INDEX_NAME)
        base_name = self.get_index_name(self.INDEX_NAME)

        replicas = [
            f"{base_name}_name_asc",
            f"{base_name}_brightness_desc",
            f"{base_name}_date_desc",
            f"{base_name}_views_desc",
        ]

        base_index.set_settings({
            'replicas': replicas
        })

        # Configure each replica
        client = self.client

        # Name A-Z
        client.init_index(f"{base_name}_name_asc").set_settings({
            'ranking': ['asc(name)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Brightness descending
        client.init_index(f"{base_name}_brightness_desc").set_settings({
            'ranking': ['desc(local_brightness)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Recently added
        client.init_index(f"{base_name}_date_desc").set_settings({
            'ranking': ['desc(created)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Most popular
        client.init_index(f"{base_name}_views_desc").set_settings({
            'ranking': ['desc(ga_views)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })


class OrganismIndexer(AlgoliaIndexer):
    """Algolia indexing for Organism model."""

    INDEX_NAME = "Organism"

    def serialize_organism(self, organism) -> dict[str, Any]:
        """Convert organism to Algolia record."""
        return {
            'objectID': str(organism.id),
            'scientific_name': organism.scientific_name,
            'division': organism.division,
            'url': organism.get_absolute_url(),
        }

    def index_organism(self, organism) -> None:
        """Index a single organism."""
        record = self.serialize_organism(organism)
        index = self.get_index(self.INDEX_NAME)
        index.save_object(record)

    def configure_index(self) -> None:
        """Configure organism index settings."""
        index = self.get_index(self.INDEX_NAME)
        index.set_settings({
            'searchableAttributes': ['scientific_name', 'division'],
        })


class ReferenceIndexer(AlgoliaIndexer):
    """Algolia indexing for Reference model."""

    INDEX_NAME = "Reference"

    def serialize_reference(self, reference) -> dict[str, Any]:
        """Convert reference to Algolia record."""
        return {
            'objectID': str(reference.id),
            'doi': reference.doi,
            'pmid': reference.pmid,
            'title': reference.title,
            'citation': reference.citation,
            'journal': reference.journal,
            'year': reference.year,
            'first_author': reference.first_author,
            'date': int(reference.date.timestamp()) if reference.date else None,
            'url': reference.get_absolute_url(),

            # Related proteins
            'prot_primary': [p.name for p in reference.primary_proteins.all()],
            'prot_secondary': [p.name for p in reference.secondary_proteins.all()],

            # Excerpts for search
            '_excerpts': reference.excerpts if hasattr(reference, 'excerpts') else [],
        }

    def index_reference(self, reference) -> None:
        """Index a single reference."""
        record = self.serialize_reference(reference)
        index = self.get_index(self.INDEX_NAME)
        index.save_object(record)

    def configure_index(self) -> None:
        """Configure reference index settings."""
        index = self.get_index(self.INDEX_NAME)
        index.set_settings({
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
        })


# Singleton instances
protein_indexer = ProteinIndexer()
organism_indexer = OrganismIndexer()
reference_indexer = ReferenceIndexer()
```

#### Step 1.3: Async Indexing with Celery

Create `backend/proteins/tasks.py`:

```python
"""Celery tasks for async Algolia indexing."""

from celery import shared_task
from django.apps import apps
import structlog

logger = structlog.get_logger(__name__)


@shared_task(bind=True, max_retries=3)
def index_protein_task(self, protein_id: int):
    """Async task to index a single protein."""
    from proteins.models import Protein
    from proteins.algolia import protein_indexer

    try:
        protein = Protein.objects.get(id=protein_id)
        protein_indexer.index_protein(protein)
        logger.info("protein_indexed", protein_id=protein_id, protein_name=protein.name)
    except Protein.DoesNotExist:
        logger.warning("protein_not_found", protein_id=protein_id)
    except Exception as exc:
        logger.error("protein_index_failed", protein_id=protein_id, error=str(exc))
        raise self.retry(exc=exc, countdown=60)  # Retry after 1 minute


@shared_task(bind=True, max_retries=3)
def delete_protein_task(self, protein_uuid: str):
    """Async task to delete a protein from index."""
    from proteins.algolia import protein_indexer
    from proteins.models import Protein

    try:
        # Create a temporary object just for deletion
        class TempProtein:
            def __init__(self, uuid):
                self.uuid = uuid

        protein_indexer.delete_protein(TempProtein(protein_uuid))
        logger.info("protein_deleted_from_index", protein_uuid=protein_uuid)
    except Exception as exc:
        logger.error("protein_delete_failed", protein_uuid=protein_uuid, error=str(exc))
        raise self.retry(exc=exc, countdown=60)


@shared_task
def reindex_all_proteins():
    """Reindex all proteins in batches."""
    from proteins.models import Protein
    from proteins.algolia import protein_indexer

    logger.info("reindex_all_proteins_started")

    proteins = Protein.objects.filter(status='approved').select_related(
        'default_state', 'primary_reference'
    )

    batch_size = 100
    total = proteins.count()

    for i in range(0, total, batch_size):
        batch = proteins[i:i + batch_size]
        protein_indexer.index_proteins_batch(batch)
        logger.info("batch_indexed", batch_start=i, batch_end=i+batch_size, total=total)

    logger.info("reindex_all_proteins_completed", total_indexed=total)


@shared_task
def configure_indices():
    """Configure all Algolia indices with proper settings."""
    from proteins.algolia import protein_indexer, organism_indexer, reference_indexer

    logger.info("configuring_algolia_indices")

    protein_indexer.configure_index()
    organism_indexer.configure_index()
    reference_indexer.configure_index()

    logger.info("algolia_indices_configured")
```

#### Step 1.4: Signal Handlers (Replace algoliasearch-django auto-indexing)

Update `backend/proteins/handlers.py`:

```python
"""Signal handlers for Algolia indexing."""

from django.db.models.signals import post_save, post_delete, m2m_changed
from django.dispatch import receiver
from django.conf import settings
from proteins.models import Protein, Organism
from references.models import Reference


# Only index if Algolia is configured
ALGOLIA_ENABLED = bool(settings.ALGOLIA.get('API_KEY'))


@receiver(post_save, sender=Protein)
def protein_saved(sender, instance, created, **kwargs):
    """Index protein when saved."""
    if not ALGOLIA_ENABLED:
        return

    # Import here to avoid circular imports
    from proteins.tasks import index_protein_task

    # Queue async indexing (doesn't block the request)
    index_protein_task.delay(instance.id)


@receiver(post_delete, sender=Protein)
def protein_deleted(sender, instance, **kwargs):
    """Remove protein from index when deleted."""
    if not ALGOLIA_ENABLED:
        return

    from proteins.tasks import delete_protein_task

    delete_protein_task.delay(str(instance.uuid))


# Similar for Organism and Reference...
@receiver(post_save, sender=Organism)
def organism_saved(sender, instance, **kwargs):
    if not ALGOLIA_ENABLED:
        return
    from proteins.algolia import organism_indexer
    organism_indexer.index_organism(instance)  # Sync for now, organisms change rarely


@receiver(post_save, sender=Reference)
def reference_saved(sender, instance, **kwargs):
    if not ALGOLIA_ENABLED:
        return
    from references.algolia import reference_indexer
    reference_indexer.index_reference(instance)
```

#### Step 1.5: Management Commands

Create `backend/proteins/management/commands/algolia_reindex.py`:

```python
"""Management command to reindex Algolia."""

from django.core.management.base import BaseCommand
from proteins.tasks import reindex_all_proteins, configure_indices


class Command(BaseCommand):
    help = 'Reindex all models in Algolia'

    def add_arguments(self, parser):
        parser.add_argument(
            '--configure',
            action='store_true',
            help='Configure index settings before reindexing',
        )

    def handle(self, *args, **options):
        if options['configure']:
            self.stdout.write('Configuring indices...')
            configure_indices.delay()
            self.stdout.write(self.style.SUCCESS('✓ Index configuration queued'))

        self.stdout.write('Reindexing proteins...')
        reindex_all_proteins.delay()
        self.stdout.write(self.style.SUCCESS('✓ Reindexing queued'))
```

#### Step 1.6: Update Settings

Update `backend/config/settings/base.py`:

```python
# REMOVE this entire section:
# if ALGOLIA["API_KEY"]:
#     INSTALLED_APPS += ["algoliasearch_django"]

# Algolia config stays the same
ALGOLIA_SUFFIX = "dev" if (DEBUG or ("staging" in env("SENTRY_PROJECT", default=""))) else "prod"
ALGOLIA_PUBLIC_KEY = env("ALGOLIA_PUBLIC_KEY", default="421b453d4f93e332ebd0c7f3ace29476")
ALGOLIA = {
    "APPLICATION_ID": "9WAWQMVNTB",
    "API_KEY": env("ALGOLIA_API_KEY", default=""),
    "INDEX_SUFFIX": ALGOLIA_SUFFIX,
}
```

---

### Phase 2: Modern Frontend (React InstantSearch)

**Goal**: Complete, beautiful search UI with autocomplete and advanced search page

#### Step 2.1: Dependencies

```json
// frontend/package.json
{
  "dependencies": {
    "@algolia/autocomplete-js": "^1.19.4",
    "react-instantsearch": "^7.17.0",
    "algoliasearch": "^5.17.0",
    "instantsearch.css": "^8.5.0"
  }
}
```

Remove:
- `"autocomplete.js": "^0.36.0"` ❌
- Downgrade `"algoliasearch": "^3.35.1"` → `"^5.17.0"` ✅

#### Step 2.2: Autocomplete Component (Navigation Bar)

Create `frontend/src/components/SearchAutocomplete.jsx`:

```jsx
import React, { createElement, Fragment, useEffect, useRef } from 'react';
import { createRoot } from 'react-dom/client';
import { autocomplete } from '@algolia/autocomplete-js';
import { createLocalStorageRecentSearchesPlugin } from '@algolia/autocomplete-plugin-recent-searches';
import algoliasearch from 'algoliasearch/lite';
import '@algolia/autocomplete-theme-classic';

const ALGOLIA_CONFIG = window.FPBASE.ALGOLIA;
const searchClient = algoliasearch(ALGOLIA_CONFIG.appID, ALGOLIA_CONFIG.publicKey);

// Recent searches plugin
const recentSearchesPlugin = createLocalStorageRecentSearchesPlugin({
  key: 'fpbase-search',
  limit: 5,
  transformSource({ source }) {
    return {
      ...source,
      templates: {
        ...source.templates,
        header() {
          return <div className="aa-SourceHeader">Recent Searches</div>;
        },
      },
    };
  },
});

function ProteinHit({ hit, components }) {
  let iconColor = 'gray50';
  if (hit.switchType && hit.switchType !== 'Basic') {
    iconColor = 'rainbow';
  } else if (hit.color && !hit.color.includes('Stokes')) {
    iconColor = hit.color.toLowerCase().replace(/ |\//g, '_');
  }

  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon aa-ItemIcon--protein">
          <img
            src={`${window.FPBASE.imageDir}gfp_${iconColor}_40.png`}
            alt=""
          />
        </div>
        <div className="aa-ItemContentBody">
          <div className="aa-ItemContentTitle">
            <components.Highlight hit={hit} attribute="name" />
          </div>
          {hit.aliases && hit.aliases.length > 0 && (
            <div className="aa-ItemContentDescription">
              aka: {hit.aliases.slice(0, 3).join(', ')}
            </div>
          )}
          {hit.ex && hit.em && (
            <div className="aa-ItemContentSubtitle">
              {hit.ex}/{hit.em} nm
            </div>
          )}
        </div>
        {hit.img_url && (
          <img
            src={hit.img_url}
            alt="spectrum"
            className="aa-ItemContentSpectra"
          />
        )}
      </div>
    </a>
  );
}

function ReferenceHit({ hit, components }) {
  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon">
          <img src={`${window.FPBASE.imageDir}ref.png`} alt="" />
        </div>
        <div className="aa-ItemContentBody">
          <div className="aa-ItemContentTitle">
            <components.Highlight hit={hit} attribute="citation" />
          </div>
          {hit._highlightResult?.title && (
            <div className="aa-ItemContentDescription">
              <components.Highlight hit={hit} attribute="title" />
            </div>
          )}
        </div>
      </div>
    </a>
  );
}

function OrganismHit({ hit, components }) {
  return (
    <a href={hit.url} className="aa-ItemLink">
      <div className="aa-ItemContent">
        <div className="aa-ItemIcon">
          <img src={`${window.FPBASE.imageDir}organism_icon.png`} alt="" />
        </div>
        <div className="aa-ItemContentBody">
          <components.Highlight hit={hit} attribute="scientific_name" />
        </div>
      </div>
    </a>
  );
}

export function SearchAutocomplete({ inputElement }) {
  const containerRef = useRef(null);
  const autocompleteRef = useRef(null);

  useEffect(() => {
    if (!containerRef.current) return;

    autocompleteRef.current = autocomplete({
      container: containerRef.current,
      placeholder: 'Search proteins, organisms, references...',
      openOnFocus: true,
      plugins: [recentSearchesPlugin],
      detachedMediaQuery: '(max-width: 768px)', // Mobile: full-screen overlay

      getSources({ query }) {
        if (!query) return [];

        return [
          // Proteins
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
                },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header() {
                return <div className="aa-SourceHeader">Proteins</div>;
              },
              item({ item, components }) {
                return <ProteinHit hit={item} components={components} />;
              },
              noResults() {
                return <div className="aa-ItemContent">No proteins found</div>;
              },
            },
            getItemUrl({ item }) {
              return item.url;
            },
          },

          // References
          {
            sourceId: 'references',
            getItems() {
              return searchClient.search([{
                indexName: ALGOLIA_CONFIG.referenceIndex,
                query,
                params: { hitsPerPage: 3 },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header() {
                return <div className="aa-SourceHeader">References</div>;
              },
              item({ item, components }) {
                return <ReferenceHit hit={item} components={components} />;
              },
            },
            getItemUrl({ item }) {
              return item.url;
            },
          },

          // Organisms
          {
            sourceId: 'organisms',
            getItems() {
              return searchClient.search([{
                indexName: ALGOLIA_CONFIG.organismIndex,
                query,
                params: { hitsPerPage: 2 },
              }]).then(({ results }) => results[0].hits);
            },
            templates: {
              header() {
                return <div className="aa-SourceHeader">Organisms</div>;
              },
              item({ item, components }) {
                return <OrganismHit hit={item} components={components} />;
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

      render({ children, state, Fragment }, root) {
        // Custom render to use React 18 createRoot
        if (!root.__root) {
          root.__root = createRoot(root);
        }
        root.__root.render(children);
      },

      renderNoResults({ state, Fragment }, root) {
        const query = state.query;
        if (!root.__root) {
          root.__root = createRoot(root);
        }
        root.__root.render(
          <div className="aa-NoResults">
            <div className="aa-NoResultsIcon">🔍</div>
            <div>No results for "{query}"</div>
            <a href={`/search/?q=${encodeURIComponent(query)}`} className="aa-NoResultsLink">
              Try advanced search →
            </a>
          </div>
        );
      },
    });

    return () => {
      autocompleteRef.current?.destroy();
    };
  }, []);

  return <div ref={containerRef} />;
}

// Export function to initialize on existing input
export function initAutocomplete() {
  const searchInput = document.querySelector('#algolia-search-input');
  if (!searchInput) {
    console.warn('Search input not found');
    return;
  }

  // Replace input with autocomplete container
  const container = document.createElement('div');
  searchInput.parentElement.replaceChild(container, searchInput);

  const root = createRoot(container);
  root.render(<SearchAutocomplete />);
}
```

#### Step 2.3: Advanced Search Page

Create `frontend/src/components/AdvancedSearch.jsx`:

```jsx
import React, { useState } from 'react';
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
  RangeInput,
  CurrentRefinements,
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
    <article className="protein-hit">
      <div className="protein-hit__icon">
        <img
          src={`${window.FPBASE.imageDir}gfp_${iconColor}_40.png`}
          alt="protein"
        />
      </div>

      <div className="protein-hit__content">
        <h3 className="protein-hit__name">
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
              λex/em: {hit.ex}/{hit.em} nm
            </span>
          )}
          {hit.switchType && hit.switchType !== 'Basic' && (
            <span className="protein-hit__badge protein-hit__badge--switch">
              {hit.switchType}
            </span>
          )}
          {hit.color && (
            <span className="protein-hit__badge protein-hit__badge--color">
              {hit.color}
            </span>
          )}
          {hit.local_brightness && (
            <span className="protein-hit__brightness">
              ε·QY: {hit.local_brightness.toFixed(1)}
            </span>
          )}
        </div>

        {hit.first_author && (
          <p className="protein-hit__author">{hit.first_author}</p>
        )}
      </div>

      {hit.img_url && (
        <div className="protein-hit__spectra">
          <img src={hit.img_url} alt="spectra" />
        </div>
      )}
    </article>
  );
}

export function AdvancedSearch() {
  const [mobileFiltersOpen, setMobileFiltersOpen] = useState(false);

  return (
    <div className="advanced-search">
      <InstantSearch
        searchClient={searchClient}
        indexName={window.FPBASE.ALGOLIA.proteinIndex}
        routing={true} // URL sync
      >
        <Configure hitsPerPage={20} />

        <div className="search-header">
          <SearchBox
            placeholder="Search proteins by name, sequence, author..."
            autoFocus
            classNames={{
              root: 'search-header__searchbox',
              input: 'search-header__input',
              submitIcon: 'search-header__submit-icon',
            }}
          />

          <div className="search-header__stats">
            <Stats
              translations={{
                rootElementText({ nbHits, processingTimeMS }) {
                  return `${nbHits.toLocaleString()} proteins found in ${processingTimeMS}ms`;
                },
              }}
            />
          </div>

          <button
            className="search-header__mobile-filter-toggle"
            onClick={() => setMobileFiltersOpen(!mobileFiltersOpen)}
          >
            {mobileFiltersOpen ? 'Hide' : 'Show'} Filters
          </button>
        </div>

        <div className="search-layout">
          <aside className={`search-sidebar ${mobileFiltersOpen ? 'search-sidebar--open' : ''}`}>
            <Panel header="Sort Results">
              <SortBy
                items={[
                  {
                    label: 'Relevance',
                    value: window.FPBASE.ALGOLIA.proteinIndex
                  },
                  {
                    label: 'Name (A-Z)',
                    value: `${window.FPBASE.ALGOLIA.proteinIndex}_name_asc`
                  },
                  {
                    label: 'Brightness',
                    value: `${window.FPBASE.ALGOLIA.proteinIndex}_brightness_desc`
                  },
                  {
                    label: 'Most Popular',
                    value: `${window.FPBASE.ALGOLIA.proteinIndex}_views_desc`
                  },
                  {
                    label: 'Recently Added',
                    value: `${window.FPBASE.ALGOLIA.proteinIndex}_date_desc`
                  },
                ]}
              />
            </Panel>

            <div className="search-sidebar__clear">
              <ClearRefinements
                translations={{
                  resetButtonText: 'Clear all filters',
                }}
              />
            </div>

            <CurrentRefinements />

            <Panel header="Switch Type" className="filter-panel">
              <RefinementList
                attribute="switchType"
                sortBy={['name:asc']}
                showMore={true}
                limit={5}
              />
            </Panel>

            <Panel header="Color" className="filter-panel">
              <RefinementList
                attribute="color"
                sortBy={['name:asc']}
                showMore={true}
                limit={10}
              />
            </Panel>

            <Panel header="Cofactor" className="filter-panel">
              <RefinementList
                attribute="cofactor"
                showMore={true}
              />
            </Panel>

            <Panel header="Brightness (ε·QY)" className="filter-panel">
              <RangeInput
                attribute="local_brightness"
                translations={{
                  separatorElementText: 'to',
                  submitButtonText: 'Apply',
                }}
              />
            </Panel>

            <Panel header="Excitation (nm)" className="filter-panel">
              <RangeInput
                attribute="ex"
                translations={{
                  separatorElementText: 'to',
                  submitButtonText: 'Apply',
                }}
              />
            </Panel>

            <Panel header="Emission (nm)" className="filter-panel">
              <RangeInput
                attribute="em"
                translations={{
                  separatorElementText: 'to',
                  submitButtonText: 'Apply',
                }}
              />
            </Panel>
          </aside>

          <main className="search-results">
            <Hits hitComponent={ProteinHit} />

            <div className="search-pagination">
              <Pagination
                showFirst={true}
                showLast={true}
                showPrevious={true}
                showNext={true}
              />
            </div>
          </main>
        </div>
      </InstantSearch>
    </div>
  );
}
```

#### Step 2.4: Entry Points

**Update `frontend/src/index.js`** (main bundle):
```javascript
import { initAutocomplete } from './components/SearchAutocomplete';

// Initialize autocomplete when DOM ready
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initAutocomplete);
} else {
  initAutocomplete();
}
```

**Create `frontend/src/search.js`** (search page bundle):
```javascript
import React from 'react';
import { createRoot } from 'react-dom/client';
import { AdvancedSearch } from './components/AdvancedSearch';

document.addEventListener('DOMContentLoaded', () => {
  const container = document.getElementById('advanced-search-root');
  if (container) {
    const root = createRoot(container);
    root.render(<AdvancedSearch />);
  }
});
```

#### Step 2.5: Update Templates

**Remove CDN script from `base.html`:**
```diff
- <script defer src="https://cdnjs.cloudflare.com/ajax/libs/autocomplete.js/0.38.1/autocomplete.jquery.min.js" ...></script>
```

**Create `backend/fpbase/templates/pages/search.html`:**
```django
{% extends "base.html" %}
{% load django_vite %}

{% block title %}Search Proteins - FPbase{% endblock %}

{% block content %}
<div id="advanced-search-root"></div>
{% vite_asset 'src/search.js' %}
{% endblock %}
```

**Update URL routing:**
```python
# backend/config/urls.py
from django.views.generic import TemplateView

urlpatterns = [
    # ... existing routes
    path('search/', TemplateView.as_view(template_name='pages/search.html'), name='search'),
]
```

---

## Benefits of This Architecture

### 1. **Frontend-First = Performance**
- **10x faster** search (direct to Algolia, no backend hop)
- Global CDN distribution (Algolia's infrastructure)
- No Django blocking on search queries
- Better mobile experience

### 2. **No Stale Dependencies**
- Direct `algoliasearch` Python client (actively maintained)
- Modern React InstantSearch (v7, latest)
- `@algolia/autocomplete-js` (actively maintained)
- No jQuery dependency for search

### 3. **Explicit Over Magic**
- You control exactly what gets indexed
- Clear serialization logic
- Easy to debug (no wrapper abstraction)
- Obvious where indexing happens

### 4. **Async by Default**
- Celery handles all indexing
- No blocking on `protein.save()`
- Batch operations for efficiency
- Resilient to Algolia API issues

### 5. **Modern UX**
- Faceted search (filters)
- Multiple sort options
- Range sliders for numeric values
- Recent searches
- Mobile-optimized
- URL routing (shareable searches)

### 6. **Maintainable**
- ~300 lines of indexing code (vs opaque wrapper)
- Standard Django patterns (signals, tasks, management commands)
- Type hints for IDE support
- Easy to extend (just add fields to `serialize_*` methods)

---

## Migration Path

### Step 1: Deploy Backend Changes (No Frontend Changes Yet)

```bash
# 1. Update dependencies
uv sync

# 2. Remove old index files
rm backend/proteins/index.py
rm backend/references/index.py

# 3. Create new files
# - backend/proteins/algolia.py
# - backend/proteins/tasks.py
# - backend/proteins/management/commands/algolia_reindex.py

# 4. Update handlers.py

# 5. Configure indices
uv run backend/manage.py shell
>>> from proteins.algolia import protein_indexer
>>> protein_indexer.configure_index()

# 6. Reindex everything
uv run backend/manage.py algolia_reindex --configure

# 7. Test indexing
uv run backend/manage.py shell
>>> from proteins.models import Protein
>>> p = Protein.objects.first()
>>> p.save()  # Should queue indexing task
```

### Step 2: Deploy Frontend Changes

```bash
# 1. Update dependencies
pnpm install

# 2. Create new components
# - frontend/src/components/SearchAutocomplete.jsx
# - frontend/src/components/AdvancedSearch.jsx

# 3. Update entry points
# - frontend/src/index.js
# - frontend/src/search.js (new)

# 4. Build
pnpm build

# 5. Test locally
pnpm dev
```

### Step 3: Gradual Rollout

1. **Week 1**: Deploy backend, keep old frontend (backwards compatible)
2. **Week 2**: Deploy new autocomplete only (test in production)
3. **Week 3**: Deploy advanced search page
4. **Week 4**: Remove old algolia.js entirely

---

## Cost Analysis

**This entire redesign costs $0**
- All open-source libraries
- No Algolia plan changes
- Same infrastructure (Celery already exists)
- Reduces bundle size (removes jQuery from search)

---

## Rollback Plan

If anything breaks:

```bash
# Backend rollback
git revert <commit>
uv sync
uv run backend/manage.py migrate

# Frontend rollback
git revert <commit>
pnpm install
pnpm build
```

Old indices remain in Algolia, so switching back is instant.

---

## Summary: Why This is Better

| Aspect | Old (algoliasearch-django) | New (Frontend-First) |
|--------|---------------------------|----------------------|
| **Backend Library** | Stale wrapper (11mo old) | Direct client (active) |
| **Frontend Library** | Deprecated autocomplete.js | Modern @algolia/autocomplete-js |
| **Search Speed** | Backend hop required | Direct to Algolia (10x faster) |
| **Indexing** | Synchronous (blocks requests) | Async via Celery |
| **Control** | Magic wrapper | Explicit serialization |
| **Code Size** | Opaque abstraction | ~300 lines, clear logic |
| **UX Features** | Basic autocomplete only | Autocomplete + faceted search |
| **Maintenance** | Depends on stale package | Direct control, easy to update |
| **React Integration** | None (jQuery) | Native React components |
| **Mobile** | Poor | Excellent (InstantSearch mobile-first) |
| **Bundle Size** | Larger (jQuery) | Smaller (tree-shakeable) |
| **Future-Proof** | No | Yes (active ecosystem) |

---

## Next Steps

**Your call:**

1. **Go all-in**: Implement the full redesign (backend + frontend)
2. **Backend first**: Just replace algoliasearch-django wrapper with custom indexing
3. **Frontend first**: Keep backend as-is, modernize UI only
4. **Hybrid**: Do backend now, frontend in phases

**My recommendation:** Start with backend (Phase 1) this week. It's low-risk, removes the stale dependency, and sets you up for the modern frontend whenever you're ready.

What do you think?
