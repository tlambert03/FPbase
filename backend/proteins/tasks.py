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
            for dye in DyeState.objects.filter(id__in=batch_ids).iterator():
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
