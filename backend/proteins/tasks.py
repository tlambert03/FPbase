from celery import shared_task
from sentry_sdk import capture_exception
import structlog

from .util.helpers import forster_list

logger = structlog.get_logger(__name__)


@shared_task
def calc_fret():
    return forster_list()


@shared_task(bind=True)
def calculate_scope_report(self, scope_id, outdated_ids=None, fluor_collection=None):
    import gc

    from proteins.models import Dye, Microscope, OcFluorEff, State

    # Initialize state_ids and dye_ids
    if not fluor_collection:
        # Use iterator to avoid loading all objects into memory at once
        # Build list of IDs instead of full objects
        state_ids = list(State.objects.with_spectra().values_list("id", flat=True))
        dye_ids = list(Dye.objects.with_spectra().values_list("id", flat=True))
    else:
        # fluor_collection is not currently implemented, but initialize to empty lists
        # to prevent potential errors if this parameter is used in the future
        state_ids = []
        dye_ids = []

    m = Microscope.objects.get(id=scope_id)
    updated = []
    i = 0
    if outdated_ids:
        total = len(outdated_ids)
        # Process in batches to control memory
        batch_size = 50
        for start in range(0, len(outdated_ids), batch_size):
            batch_ids = outdated_ids[start : start + batch_size]
            for x in OcFluorEff.objects.filter(id__in=batch_ids):
                i += 1
                self.update_state(state="PROGRESS", meta={"current": i, "total": total})
                x.save()
            # Force garbage collection after each batch
            gc.collect()
        return

    # Process states and dyes separately in batches to reduce memory usage
    oc_count = m.optical_configs.count()
    total = oc_count * (len(state_ids) + len(dye_ids))

    # Process states in batches
    batch_size = 50
    for oc in m.optical_configs.all():
        # Process State objects in batches
        for start in range(0, len(state_ids), batch_size):
            batch_ids = state_ids[start : start + batch_size]
            for state in State.objects.filter(id__in=batch_ids).iterator():
                i += 1
                self.update_state(state="PROGRESS", meta={"current": i, "total": total})
                try:
                    obj = OcFluorEff.objects.get(oc=oc, state=state)
                    if obj.outdated:
                        obj.save()
                        updated.append((oc, state))
                except OcFluorEff.DoesNotExist:
                    try:
                        OcFluorEff.objects.create(oc=oc, fluor=state)
                    except Exception as e:
                        capture_exception(e)
            gc.collect()

        # Process Dye objects in batches
        for start in range(0, len(dye_ids), batch_size):
            batch_ids = dye_ids[start : start + batch_size]
            for dye in Dye.objects.filter(id__in=batch_ids).iterator():
                i += 1
                self.update_state(state="PROGRESS", meta={"current": i, "total": total})
                try:
                    obj = OcFluorEff.objects.get(oc=oc, dye=dye)
                    if obj.outdated:
                        obj.save()
                        updated.append((oc, dye))
                except OcFluorEff.DoesNotExist:
                    try:
                        OcFluorEff.objects.create(oc=oc, fluor=dye)
                    except Exception as e:
                        capture_exception(e)
            gc.collect()


# ============================================================================
# Algolia Indexing Tasks
# ============================================================================


@shared_task(bind=True, max_retries=3)
def index_protein_task(self, protein_id: int):
    """Async task to index a single protein in Algolia."""
    from proteins.models import Protein
    from proteins.algolia import protein_indexer

    try:
        protein = Protein.objects.select_related(
            'default_state', 'primary_reference'
        ).get(id=protein_id)
        protein_indexer.index_protein(protein)
        logger.info("protein_indexing_task_completed", protein_id=protein_id, protein_name=protein.name)
    except Protein.DoesNotExist:
        logger.warning("protein_not_found_for_indexing", protein_id=protein_id)
    except Exception as exc:
        logger.error("protein_indexing_task_failed", protein_id=protein_id, error=str(exc), exc_info=True)
        # Retry with exponential backoff (60s, 120s, 240s)
        raise self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))


@shared_task(bind=True, max_retries=3)
def delete_protein_task(self, protein_uuid: str):
    """Async task to delete a protein from Algolia index."""
    from proteins.algolia import protein_indexer

    try:
        # Create a temporary object just for deletion
        class TempProtein:
            def __init__(self, uuid):
                self.uuid = uuid

        protein_indexer.delete_protein(TempProtein(protein_uuid))
        logger.info("protein_deletion_task_completed", protein_uuid=protein_uuid)
    except Exception as exc:
        logger.error("protein_deletion_task_failed", protein_uuid=protein_uuid, error=str(exc))
        raise self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))


@shared_task
def reindex_all_proteins():
    """Reindex all approved proteins in batches."""
    from proteins.models import Protein
    from proteins.algolia import protein_indexer

    logger.info("reindex_all_proteins_started")

    try:
        proteins = Protein.objects.filter(
            status='approved'
        ).select_related(
            'default_state', 'primary_reference'
        ).order_by('id')

        batch_size = 100
        total = proteins.count()

        logger.info("reindex_total_count", total=total)

        for i in range(0, total, batch_size):
            batch = list(proteins[i:i + batch_size])
            protein_indexer.index_proteins_batch(batch)
            logger.info("batch_indexed", batch_start=i, batch_size=len(batch), total=total)

        logger.info("reindex_all_proteins_completed", total_indexed=total)
        return {'total_indexed': total}

    except Exception as exc:
        logger.error("reindex_all_proteins_failed", error=str(exc), exc_info=True)
        raise


@shared_task
def configure_protein_index():
    """Configure protein index with proper settings and replicas."""
    from proteins.algolia import protein_indexer

    logger.info("configuring_protein_index")

    try:
        protein_indexer.configure_index()
        logger.info("protein_index_configured")
        return {'status': 'configured'}
    except Exception as exc:
        logger.error("protein_index_configuration_failed", error=str(exc), exc_info=True)
        raise


@shared_task
def index_organism_task(organism_id: int):
    """Async task to index a single organism in Algolia."""
    from proteins.models import Organism
    from proteins.algolia import organism_indexer

    try:
        organism = Organism.objects.get(id=organism_id)
        organism_indexer.index_organism(organism)
        logger.info("organism_indexed", organism_id=organism_id)
    except Organism.DoesNotExist:
        logger.warning("organism_not_found", organism_id=organism_id)
    except Exception as exc:
        logger.error("organism_indexing_failed", organism_id=organism_id, error=str(exc))
        raise


@shared_task
def index_reference_task(reference_id: int):
    """Async task to index a single reference in Algolia."""
    from references.models import Reference
    from proteins.algolia import reference_indexer

    try:
        reference = Reference.objects.prefetch_related(
            'primary_proteins', 'proteins'
        ).get(id=reference_id)
        reference_indexer.index_reference(reference)
        logger.info("reference_indexed", reference_id=reference_id)
    except Reference.DoesNotExist:
        logger.warning("reference_not_found", reference_id=reference_id)
    except Exception as exc:
        logger.error("reference_indexing_failed", reference_id=reference_id, error=str(exc))
        raise


@shared_task
def configure_all_indices():
    """Configure all Algolia indices with proper settings."""
    from proteins.algolia import protein_indexer, organism_indexer, reference_indexer

    logger.info("configuring_all_indices")

    results = {}

    try:
        protein_indexer.configure_index()
        results['protein'] = 'configured'
    except Exception as exc:
        logger.error("protein_index_config_failed", error=str(exc))
        results['protein'] = f'failed: {exc}'

    try:
        organism_indexer.configure_index()
        results['organism'] = 'configured'
    except Exception as exc:
        logger.error("organism_index_config_failed", error=str(exc))
        results['organism'] = f'failed: {exc}'

    try:
        reference_indexer.configure_index()
        results['reference'] = 'configured'
    except Exception as exc:
        logger.error("reference_index_config_failed", error=str(exc))
        results['reference'] = f'failed: {exc}'

    logger.info("all_indices_configured", results=results)
    return results
