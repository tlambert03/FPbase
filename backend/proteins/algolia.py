"""
Custom Algolia indexing service for FPbase.

Why not algoliasearch-django?
- That package is inactive/discontinued (last release 7mo ago)
- Frontend-first architecture means minimal backend indexing needs
- Direct client gives us more control and is actively maintained
- Explicit serialization is easier to understand and debug

This module provides custom indexing services for Protein, Organism, and Reference models.
All indexing is done asynchronously via Celery tasks to avoid blocking request/response cycles.
"""

from typing import Any
from django.conf import settings
from algoliasearch.search.client import SearchClient
import structlog

logger = structlog.get_logger(__name__)


class AlgoliaIndexer:
    """Base indexing service for Algolia."""

    def __init__(self):
        if not settings.ALGOLIA.get('API_KEY'):
            logger.warning("algolia_api_key_not_configured")
            self.client = None
            return

        self.client = SearchClient.create(
            settings.ALGOLIA['APPLICATION_ID'],
            settings.ALGOLIA['API_KEY']
        )
        self.index_suffix = settings.ALGOLIA['INDEX_SUFFIX']

    def get_index_name(self, base_name: str) -> str:
        """Get full index name with environment suffix (dev/prod)."""
        return f"{base_name}_{self.index_suffix}"

    def get_index(self, base_name: str):
        """Get index instance."""
        if not self.client:
            return None
        return self.client.init_index(self.get_index_name(base_name))


class ProteinIndexer(AlgoliaIndexer):
    """Algolia indexing for Protein model."""

    INDEX_NAME = "Protein"

    def serialize_protein(self, protein) -> dict[str, Any] | None:
        """Convert protein to Algolia record."""
        if not protein.is_visible():
            return None  # Don't index hidden proteins

        # Get the display values for choice fields
        try:
            switch_type_display = protein.get_switch_type_display() if protein.switch_type else None
        except Exception:
            switch_type_display = None

        try:
            agg_display = protein.get_agg_display() if protein.agg else None
        except Exception:
            agg_display = None

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

            # External IDs
            'pdb': list(protein.pdb) if protein.pdb else [],
            'genbank': protein.genbank,
            'uniprot': protein.uniprot,
            'ipg_id': protein.ipg_id,

            # Properties for faceting
            'switchType': switch_type_display,
            'agg': agg_display,

            # Reference info
            'first_author': protein.first_author() if hasattr(protein, 'first_author') else None,

            # Popularity metrics for custom ranking
            'rank': getattr(protein, 'rank', 0) or 0,
        }

        # Add default_state properties if available
        if protein.default_state:
            record.update({
                'ex': protein.default_state.ex_max,
                'em': protein.default_state.em_max,
                'pka': protein.default_state.pka,
                'ec': protein.default_state.ext_coeff,
                'qy': protein.default_state.qy,
            })

            # Get color from color property
            try:
                record['color'] = protein.color() if hasattr(protein, 'color') else None
            except Exception:
                record['color'] = None

            # Get em_css if available
            try:
                record['em_css'] = protein.em_css() if hasattr(protein, 'em_css') else None
            except Exception:
                record['em_css'] = None

        # Add cofactor
        try:
            record['cofactor'] = protein.get_cofactor_display() if protein.cofactor else None
        except Exception:
            record['cofactor'] = None

        # Computed brightness
        try:
            record['local_brightness'] = protein.local_brightness() if hasattr(protein, 'local_brightness') else None
        except Exception:
            record['local_brightness'] = None

        # Analytics/popularity data - these might not exist on all instances
        record['ga_views'] = getattr(protein, 'ga_views', 0) or 0
        record['n_faves'] = getattr(protein, 'n_faves', 0) or 0
        record['n_cols'] = getattr(protein, 'n_cols', 0) or 0

        # Dates
        if protein.created:
            record['created'] = int(protein.created.timestamp())

        if protein.primary_reference and protein.primary_reference.date:
            record['date_published'] = int(protein.primary_reference.date.timestamp())

        # UI helper - spectrum image URL
        try:
            if hasattr(protein, 'img_url'):
                img = protein.img_url()
                record['img_url'] = img if img else None
        except Exception:
            record['img_url'] = None

        # Tags for filtering
        record['_tags'] = self._get_tags(protein)

        # Remove None values to keep index clean
        return {k: v for k, v in record.items() if v is not None}

    def _get_tags(self, protein) -> list[str]:
        """Get tags for filtering."""
        tags = []

        if protein.switch_type and protein.switch_type != 'b':  # Not basic
            tags.append('switchable')

        if protein.agg == 'm':
            tags.append('monomer')

        try:
            color = protein.color() if hasattr(protein, 'color') else None
            if color:
                tags.append(f'color:{color}')
        except Exception:
            pass

        return tags

    def index_protein(self, protein) -> None:
        """Index a single protein."""
        if not self.client:
            logger.debug("algolia_not_configured", action="skip_index")
            return

        try:
            record = self.serialize_protein(protein)
            if record:
                index = self.get_index(self.INDEX_NAME)
                index.save_object(record)
                logger.info("protein_indexed", protein_id=protein.id, protein_name=protein.name)
            else:
                # Protein shouldn't be indexed (hidden), delete if exists
                self.delete_protein(protein)
        except Exception as exc:
            logger.error("protein_index_failed", protein_id=protein.id, error=str(exc), exc_info=True)
            raise

    def index_proteins_batch(self, proteins) -> None:
        """Index multiple proteins in batch."""
        if not self.client:
            logger.debug("algolia_not_configured", action="skip_batch_index")
            return

        try:
            records = []
            for protein in proteins:
                record = self.serialize_protein(protein)
                if record:
                    records.append(record)

            if records:
                index = self.get_index(self.INDEX_NAME)
                index.save_objects(records)
                logger.info("proteins_batch_indexed", count=len(records))
        except Exception as exc:
            logger.error("protein_batch_index_failed", error=str(exc), exc_info=True)
            raise

    def delete_protein(self, protein) -> None:
        """Remove protein from index."""
        if not self.client:
            return

        try:
            index = self.get_index(self.INDEX_NAME)
            index.delete_object(str(protein.uuid))
            logger.info("protein_deleted_from_index", protein_uuid=str(protein.uuid))
        except Exception as exc:
            logger.error("protein_delete_failed", protein_uuid=str(protein.uuid), error=str(exc))
            raise

    def clear_index(self) -> None:
        """Clear all records from index."""
        if not self.client:
            return

        index = self.get_index(self.INDEX_NAME)
        index.clear_objects()
        logger.info("protein_index_cleared")

    def configure_index(self) -> None:
        """Configure index settings (searchable attributes, ranking, etc)."""
        if not self.client:
            logger.warning("algolia_not_configured", action="skip_configure")
            return

        index = self.get_index(self.INDEX_NAME)

        settings = {
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
                'desc(ga_views)',         # Most viewed first
                'desc(n_faves)',          # Most favorited
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
        }

        index.set_settings(settings)
        logger.info("protein_index_configured", index_name=self.get_index_name(self.INDEX_NAME))

        # Create replica indices for different sort orders
        self._create_replicas()

    def _create_replicas(self) -> None:
        """Create replica indices for alternative sorting."""
        if not self.client:
            return

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
        # Name A-Z
        self.client.init_index(f"{base_name}_name_asc").set_settings({
            'ranking': ['asc(name)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Brightness descending
        self.client.init_index(f"{base_name}_brightness_desc").set_settings({
            'ranking': ['desc(local_brightness)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Recently added
        self.client.init_index(f"{base_name}_date_desc").set_settings({
            'ranking': ['desc(created)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        # Most popular
        self.client.init_index(f"{base_name}_views_desc").set_settings({
            'ranking': ['desc(ga_views)', 'typo', 'geo', 'words', 'filters', 'proximity', 'attribute', 'exact']
        })

        logger.info("protein_replicas_created", replicas=replicas)


class OrganismIndexer(AlgoliaIndexer):
    """Algolia indexing for Organism model."""

    INDEX_NAME = "Organism"

    def serialize_organism(self, organism) -> dict[str, Any]:
        """Convert organism to Algolia record."""
        return {
            'objectID': str(organism.id),
            'scientific_name': organism.scientific_name,
            'division': organism.division if hasattr(organism, 'division') else None,
            'url': organism.get_absolute_url(),
        }

    def index_organism(self, organism) -> None:
        """Index a single organism."""
        if not self.client:
            return

        try:
            record = self.serialize_organism(organism)
            index = self.get_index(self.INDEX_NAME)
            index.save_object(record)
            logger.info("organism_indexed", organism_id=organism.id)
        except Exception as exc:
            logger.error("organism_index_failed", organism_id=organism.id, error=str(exc))
            raise

    def configure_index(self) -> None:
        """Configure organism index settings."""
        if not self.client:
            return

        index = self.get_index(self.INDEX_NAME)
        index.set_settings({
            'searchableAttributes': ['scientific_name', 'division'],
        })
        logger.info("organism_index_configured")


class ReferenceIndexer(AlgoliaIndexer):
    """Algolia indexing for Reference model."""

    INDEX_NAME = "Reference"

    def serialize_reference(self, reference) -> dict[str, Any]:
        """Convert reference to Algolia record."""
        record = {
            'objectID': str(reference.id),
            'doi': reference.doi if hasattr(reference, 'doi') else None,
            'pmid': reference.pmid if hasattr(reference, 'pmid') else None,
            'title': reference.title if hasattr(reference, 'title') else None,
            'citation': reference.citation if hasattr(reference, 'citation') else None,
            'journal': reference.journal if hasattr(reference, 'journal') else None,
            'year': reference.year if hasattr(reference, 'year') else None,
            'url': reference.get_absolute_url(),
        }

        # Add first author if available
        if hasattr(reference, 'first_author') and reference.first_author:
            if hasattr(reference.first_author, 'family'):
                record['first_author'] = reference.first_author.family
            else:
                record['first_author'] = str(reference.first_author)

        # Add date timestamp if available
        if hasattr(reference, 'date') and reference.date:
            record['date'] = int(reference.date.timestamp())

        # Related proteins
        if hasattr(reference, 'primary_proteins'):
            record['prot_primary'] = [p.name for p in reference.primary_proteins.all()]

        if hasattr(reference, 'proteins'):
            record['prot_secondary'] = [p.name for p in reference.proteins.all()]

        # Excerpts for search (if available)
        if hasattr(reference, 'excerpts'):
            record['_excerpts'] = reference.excerpts

        # Remove None values
        return {k: v for k, v in record.items() if v is not None}

    def index_reference(self, reference) -> None:
        """Index a single reference."""
        if not self.client:
            return

        try:
            record = self.serialize_reference(reference)
            index = self.get_index(self.INDEX_NAME)
            index.save_object(record)
            logger.info("reference_indexed", reference_id=reference.id)
        except Exception as exc:
            logger.error("reference_index_failed", reference_id=reference.id, error=str(exc))
            raise

    def configure_index(self) -> None:
        """Configure reference index settings."""
        if not self.client:
            return

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
        logger.info("reference_index_configured")


# Singleton instances
protein_indexer = ProteinIndexer()
organism_indexer = OrganismIndexer()
reference_indexer = ReferenceIndexer()
