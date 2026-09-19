# ruff: noqa: F821  (PARAMS and emit are injected by remote.py / _bootstrap.py)
# READ-ONLY.  Deterministic triage of every pending protein: what is the NET change
# between the record as it was before it went pending and the record today?
# No literature, no judgement -- just "what did the submission actually touch".
import json
from datetime import timedelta

from django.core import serializers
from reversion.models import Revision, Version

from proteins.models import BleachMeasurement, Protein, State, StateTransition
from references.models import Reference

PROTEIN_FIELDS = (
    "name", "aliases", "seq", "seq_comment", "pdb", "genbank", "uniprot", "ipg_id", "mw",
    "agg", "oser", "switch_type", "blurb", "cofactor", "chromophore", "parent_organism",
    "primary_reference", "references",
)  # fmt: skip
STATE_FIELDS = (
    "name", "ex_max", "em_max", "ext_coeff", "qy", "pka", "lifetime", "maturation",
    "twop_ex_max", "twop_peak_gm", "twop_qy", "is_dark",
)  # fmt: skip
LINEAGE_FIELDS = ("parent", "mutation")
TRANSITION_FIELDS = ("trans_wave", "from_state", "to_state")
RELATED = ("transitions", "excerpts", "oser_measurements")  # reverse FKs on Protein
FLUOR_STATE = State._meta.get_parent_list()[0]


def current(obj):
    return json.loads(serializers.serialize("json", [obj]))[0]["fields"]


def norm(v):
    # "", None, [] and "[]" all mean "empty"; floats stored as strings in old snapshots
    if v in ("", None, [], "[]", "{}", {}):
        return None
    if isinstance(v, list):
        return sorted(map(str, v))
    try:
        return round(float(v), 4)
    except (TypeError, ValueError):
        return str(v).strip()


def snapshot_before(model, pk, cutoff):
    v = (
        Version.objects.get_for_object_reference(model, pk)
        .filter(revision__date_created__lt=cutoff)
        .order_by("-pk")
        .first()
    )
    return json.loads(v.serialized_data)[0]["fields"] if v else None


def diff(before, after, keys):
    # keys missing from an old snapshot (schema drift) are not evidence of a change
    return {
        k: [before[k], after.get(k)]
        for k in keys
        if k in before and norm(before[k]) != norm(after.get(k))
    }


def deleted_since(model, p, cutoff, current_ids):
    """Objects of `model` that belonged to `p` before `cutoff` and no longer exist."""
    versions = Version.objects.get_for_model(model).filter(
        revision__date_created__lt=cutoff, serialized_data__contains=f'"protein": {p.pk}'
    )
    gone = {}
    for v in versions.order_by("pk"):
        # `contains` is only a prefilter: '"protein": 6' also matches '"protein": 61'
        if json.loads(v.serialized_data)[0]["fields"].get("protein") != p.pk:
            continue
        if int(v.object_id) not in current_ids:
            gone[v.object_id] = v.object_repr
    return sorted(gone.values())


def ref_info(ref):
    return {k: getattr(ref, k) for k in ("id", "doi", "citation", "title", "year")}


def editors_since(p, cutoff):
    revs = Revision.objects.filter(
        id__in=p.versions.values_list("revision_id", flat=True), date_created__gte=cutoff
    ).select_related("user")
    return sorted(
        {f"{r.user.username}{' (staff)' if r.user.is_staff else ''}" for r in revs if r.user}
    )


def lineage_matches_seq(p):
    lin = getattr(p, "lineage", None)
    if not (lin and lin.parent and lin.mutation and p.seq and lin.parent.protein.seq):
        return None
    try:
        return str(lin.parent.protein.seq.mutate(lin.mutation)) == str(p.seq)
    except Exception as e:
        return f"error: {type(e).__name__}"


def triage(p):
    never_approved = (p.status_changed - p.created) < timedelta(minutes=1)
    cutoff = p.created if never_approved else p.status_changed - timedelta(minutes=1)
    out = {
        "slug": p.slug,
        "name": p.name,
        "modified": p.modified,
        "pending_since": p.status_changed,
        "last_editor": p.updated_by and p.updated_by.username,
        "last_editor_staff": bool(p.updated_by and p.updated_by.is_staff),
        "lineage_matches_seq": lineage_matches_seq(p),
        "n_states": p.states.count(),
        "has_seq": bool(p.seq),
        "primary_doi": p.primary_reference and p.primary_reference.doi,
        "editors": editors_since(p, cutoff),
        # every paper the record cites: a submitted value may come from any of them
        "cited": [ref_info(r) for r in p.references.all()],
    }
    if never_approved and p.last_approved_version() is None:
        return {**out, "kind": "new"}

    base = snapshot_before(Protein, p.pk, cutoff)
    if base is None:
        return {**out, "kind": "edit", "baseline": None}

    changes = {f"protein.{k}": v for k, v in diff(base, current(p), PROTEIN_FIELDS).items()}
    if "protein.references" in changes:  # ids -> which papers were added / removed
        old, new = (set(x or []) for x in changes.pop("protein.references"))
        refs = {r.id: r for r in Reference.objects.filter(id__in=old ^ new)}
        if new - old:
            changes["references.added"] = [
                None,
                [ref_info(refs[i]) for i in new - old if i in refs],
            ]
        if old - new:
            changes["references.removed"] = [
                [ref_info(refs[i]) for i in old - new if i in refs],
                None,
            ]
    for s in p.states.all():
        parent = s.fluorstate_ptr
        # pre-refactor snapshots keep the spectral fields on State, newer ones on FluorState
        b = {
            **(snapshot_before(State, s.pk, cutoff) or {}),
            **(snapshot_before(FLUOR_STATE, s.pk, cutoff) or {}),
        }
        if not b:
            if s.created >= cutoff:
                changes[f"state[{s.name}]"] = ["<added>", {k: getattr(s, k) for k in STATE_FIELDS}]
            continue  # never snapshotted: can't tell, don't invent a change
        now = {**current(s), **current(parent)}
        changes.update({f"state[{s.name}].{k}": v for k, v in diff(b, now, STATE_FIELDS).items()})

    lin = getattr(p, "lineage", None)
    if lin:
        b = snapshot_before(type(lin), lin.pk, cutoff)
        if b is None and lin.created >= cutoff:
            changes["lineage"] = ["<added>", f"{lin.parent} + {lin.mutation}"]
        elif b:
            changes.update(
                {f"lineage.{k}": v for k, v in diff(b, current(lin), LINEAGE_FIELDS).items()}
            )

    for t in p.transitions.all():
        b = snapshot_before(StateTransition, t.pk, cutoff)
        if b:
            changes.update(
                {
                    f"transition[{t}].{k}": v
                    for k, v in diff(b, current(t), TRANSITION_FIELDS).items()
                }
            )
    for model, ids, label in (
        (State, set(p.states.values_list("id", flat=True)), "states"),
        (StateTransition, set(p.transitions.values_list("id", flat=True)), "transitions"),
    ):
        if gone := deleted_since(model, p, cutoff, ids):
            changes[f"{label}.deleted"] = [gone, None]

    for rel in RELATED:
        n = getattr(p, rel).filter(created__gte=cutoff).count()
        if n:
            changes[rel] = ["<added>", n]
    n = BleachMeasurement.objects.filter(state__protein=p, created__gte=cutoff).count()
    if n:
        changes["bleach_measurements"] = ["<added>", n]
    return {**out, "kind": "edit", "baseline": True, "changes": changes}


qs = Protein.objects.filter(status="pending").select_related(
    "updated_by", "primary_reference", "lineage__parent__protein"
)
if PARAMS.get("slugs"):
    qs = qs.filter(slug__in=PARAMS["slugs"])

rows = []
for p in qs:
    try:
        rows.append(triage(p))
    except Exception as e:
        rows.append({"slug": p.slug, "kind": "error", "error": f"{type(e).__name__}: {e}"})
emit(rows)
