from celery import shared_task
from sentry_sdk import capture_exception

from .util.helpers import forster_list


@shared_task
def calc_fret():
    return forster_list()


@shared_task(bind=True)
def calculate_scope_report(self, scope_id, outdated_ids=None, fluor_collection=None):
    import gc

    from proteins.models import DyeState, Microscope, OcFluorEff, State

    # Initialize state_ids and dye_ids
    if not fluor_collection:
        # Use iterator to avoid loading all objects into memory at once
        # Build list of IDs instead of full objects
        state_ids = list(State.objects.with_spectra().values_list("id", flat=True))
        dye_ids = list(DyeState.objects.with_spectra().values_list("id", flat=True))
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
    # Note: We query concrete subclasses (State, DyeState) instead of Fluorophore
    # because MTI requires this to access child-specific fields
    oc_count = m.optical_configs.count()
    total = oc_count * (len(state_ids) + len(dye_ids))

    batch_size = 50
    for oc in m.optical_configs.all():
        # Process both State and DyeState with the same logic
        for model_class, fluor_ids in [(State, state_ids), (DyeState, dye_ids)]:
            for start in range(0, len(fluor_ids), batch_size):
                batch_ids = fluor_ids[start : start + batch_size]
                for fluor in model_class.objects.filter(id__in=batch_ids).iterator():
                    i += 1
                    self.update_state(state="PROGRESS", meta={"current": i, "total": total})
                    try:
                        obj = OcFluorEff.objects.get(oc=oc, fluor=fluor)
                        if obj.outdated:
                            obj.save()
                            updated.append((oc, fluor))
                    except OcFluorEff.DoesNotExist:
                        try:
                            OcFluorEff.objects.create(oc=oc, fluor=fluor)
                        except Exception as e:
                            capture_exception(e)
                gc.collect()
