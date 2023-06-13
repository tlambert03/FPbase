from celery import shared_task
from sentry_sdk import capture_exception

from .util.helpers import forster_list


@shared_task
def calc_fret():
    return forster_list()


@shared_task(bind=True)
def calculate_scope_report(self, scope_id, outdated_ids=None, fluor_collection=None):
    from proteins.models import Dye, Microscope, OcFluorEff, State

    if not fluor_collection:
        fluor_collection = list(State.objects.with_spectra())
        fluor_collection += list(Dye.objects.with_spectra())
    m = Microscope.objects.get(id=scope_id)
    updated = []
    i = 0
    if outdated_ids:
        total = len(outdated_ids)
        for x in OcFluorEff.objects.filter(id__in=outdated_ids).all():
            i += 1
            self.update_state(state="PROGRESS", meta={"current": i, "total": total})
            x.save()
        return
    total = m.optical_configs.count() * len(fluor_collection)
    for oc in m.optical_configs.all():
        for f in set(fluor_collection):
            i += 1
            self.update_state(state="PROGRESS", meta={"current": i, "total": total})
            try:
                kwargs = {"oc": oc}
                kwargs[f.__class__.__name__.lower()] = f
                obj = OcFluorEff.objects.get(**kwargs)
                if obj.outdated:
                    obj.save()
                    updated.append((oc, f))
            except OcFluorEff.DoesNotExist:
                try:
                    OcFluorEff.objects.create(oc=oc, fluor=f)
                except Exception as e:
                    capture_exception(e)
