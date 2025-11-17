# myapp/handlers.py
from corsheaders.signals import check_request_enabled
from django.conf import settings
from django.db.models.signals import post_delete, post_save
from django.dispatch import receiver


def cors_allow_api_to_everyone(sender, request, **kwargs):
    return request.path.startswith("/api/")


check_request_enabled.connect(cors_allow_api_to_everyone)


# ============================================================================
# Algolia Indexing Signal Handlers
# ============================================================================

# Only set up Algolia signals if API key is configured
ALGOLIA_ENABLED = bool(settings.ALGOLIA.get("API_KEY"))


@receiver(post_save, sender="proteins.Protein")
def protein_saved(sender, instance, created, **kwargs):
    """Index protein in Algolia when saved."""
    if not ALGOLIA_ENABLED:
        return

    # Import here to avoid circular imports
    from proteins.tasks import index_protein_task

    # Queue async indexing (doesn't block the request)
    index_protein_task.delay(instance.id)


@receiver(post_delete, sender="proteins.Protein")
def protein_deleted(sender, instance, **kwargs):
    """Remove protein from Algolia index when deleted."""
    if not ALGOLIA_ENABLED:
        return

    from proteins.tasks import delete_protein_task

    # Queue async deletion
    delete_protein_task.delay(str(instance.uuid))


@receiver(post_save, sender="proteins.Organism")
def organism_saved(sender, instance, **kwargs):
    """Index organism in Algolia when saved."""
    if not ALGOLIA_ENABLED:
        return

    from proteins.tasks import index_organism_task

    # Organisms change rarely, so we can afford to do this async
    index_organism_task.delay(instance.id)


@receiver(post_save, sender="references.Reference")
def reference_saved(sender, instance, **kwargs):
    """Index reference in Algolia when saved."""
    if not ALGOLIA_ENABLED:
        return

    from proteins.tasks import index_reference_task

    # Queue async indexing
    index_reference_task.delay(instance.id)
