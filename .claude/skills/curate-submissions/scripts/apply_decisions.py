# ruff: noqa: F821  (PARAMS and emit are injected by remote.py / _bootstrap.py)
# WRITES TO PRODUCTION when PARAMS["commit"] is true; otherwise everything is rolled back.
# Runs remotely (see remote.py).  Mirrors views.ajax.approve_protein,
# views.protein.revert_revision and views.spectra.pending_spectrum_action.
#
# PARAMS: commit (bool), moderator (staff username), decisions: list of
#   {"kind": "protein", "slug": ..., "action": "approve"|"reject"|"correct",
#    # "correct": the same edits on a record of any status, without changing its status
#    "expect_modified": <"modified" from fetch/triage>, "reason": ...,
#    # optional moderator corrections, approve only, saved in the approval's revision.
#    # "fix+approve" sets verified values; "undo edit" puts the pre-pending values back.
#    # fluorescence values live on measurements (what one paper reported); the state is
#    # the composite shown on the page.  Put each value on the row of the paper that
#    # reported it; "doi": null = source unknown.  A value of null clears the field; a
#    # row left empty is deleted.  The state is rebuilt from its measurements.
#    "measurements": [{"state": <name>, "doi": <doi or null>, "values": {<field>: <value>},
#                      "conditions": <optional text>}],
#    "state_edits": {<state name>: {<field>: <value>}},  # only name / maturation
#    "protein_edits": {<field>: <value>},
#    "lineage_mutation": "K69E/C134W/M205I",  # only accepted if parent + it == seq
#    "primary_doi": "10.1126/...",  # e.g. preprint -> published; old one stays as a reference
#    "remove_references": [<doi>, ...]}
#   {"kind": "spectrum", "id": ..., "action": "approve"|"reject", "reason": ...}
import contextlib
from datetime import timedelta

import reversion
from django.contrib.auth import get_user_model
from django.db import transaction
from django.http import HttpRequest
from reversion.models import Revision

from fpbase.util import uncache_protein_page
from proteins.models import Lineage, Protein, Spectrum, State
from references.models import Reference

TAG = "[curate-submissions]"
# the only fields a decision may edit (never status, slug, ownership, ...)
PROTEIN_EDITABLE = {
    "name", "aliases", "seq", "seq_validated", "seq_comment", "pdb", "genbank", "uniprot",
    "ipg_id", "mw", "agg", "oser", "switch_type", "blurb", "cofactor", "chromophore",
    "parent_organism_id",
}  # fmt: skip
STATE_EDITABLE = {"name", "maturation"}  # not measurements
MEASURED = set(State._written_through_fields())


class Skip(Exception):
    pass


def snapshot(p):
    def fields(obj):
        return {f.attname: f.value_from_object(obj) for f in obj._meta.concrete_fields}

    skip = {"modified", "status_changed", "updated_by_id"}
    snap = {f"protein.{k}": v for k, v in fields(p).items() if k not in skip}
    for s in p.states.all():
        snap.update({f"state[{s.name}].{k}": v for k, v in fields(s).items() if k not in skip})
        for m in s.measurements.select_related("reference"):
            doi = m.reference.doi if m.reference else "no reference"
            snap.update(
                {f"measurement[{s.name} | {doi}].{k}": getattr(m, k) for k in sorted(MEASURED)}
            )
    snap["references"] = sorted(p.references.values_list("doi", flat=True))
    snap["transitions"] = sorted(str(t) for t in p.transitions.all())
    if hasattr(p, "lineage"):
        snap["lineage"] = f"{p.lineage.parent} {p.lineage.mutation}"
    return snap


def pending_fields(p, approved):
    """Names of fields touched by revisions made since the last approved one."""
    names = {"status"}
    for rev, changes in p.history().items():
        if rev.id > approved.revision_id:
            for items in changes.values():
                # field is None when a whole related object was added/removed
                names.update(field or "<object>" for _, field, _ in items)
    return names


def collateral(diff, allowed):
    """Keys in a revert's diff that the pending revisions don't explain.

    Old revisions predate schema changes, so reverting to them can clobber unrelated data.
    """
    bad = []
    for key, (before, after) in diff.items():
        name = key.rsplit(".", 1)[-1]
        if name == "created":  # reversion truncates microseconds
            continue
        if name in ("lineage", "transitions") and "<object>" in allowed:
            continue
        retyped = None not in (before, after) and type(before) is not type(after)
        if name not in allowed or retyped:
            bad.append(key)
    return bad


def set_fields(obj, values, editable):
    for field, value in values.items():
        if field in MEASURED and field not in editable:
            raise Skip(f"{field!r} is a measurement: use `measurements`, not a state edit")
        if field not in editable:
            raise Skip(f"{field!r} is not an editable field")
        setattr(obj, field, value)


def set_measurement(p, user, m):
    """Create/update the measurement for (state, reference); the state rebuilds itself."""
    state = p.states.get(name=m["state"])
    # same call the public form makes for a DOI (looks it up on first use)
    ref = Reference.objects.get_or_create(doi=m["doi"].lower())[0] if m.get("doi") else None
    meas = state.measurements.filter(reference=ref).first()
    if meas is None:
        meas = state.measurements.model(state=state, reference=ref, created_by=user)
    set_fields(meas, m["values"], MEASURED)
    if "conditions" in m:
        meas.conditions = m["conditions"]
    meas.updated_by = user
    if any(getattr(meas, f) not in (None, False) for f in MEASURED):
        meas.save()
    elif meas.pk:
        meas.delete()


def approve_protein(p, user, reason, d, approve=True):
    if approve:
        # get rid of previous unapproved version (as in views.ajax.approve_protein)
        with contextlib.suppress(Exception):
            if p.versions.first().field_dict["status"] == "pending":
                p.versions.first().delete()
    with reversion.create_revision():
        reversion.set_user(user)
        what = "approved current version" if approve else "corrected"
        reversion.set_comment(f"{user} {what} {TAG} {reason}")
        # moderator corrections, saved in the same revision as the approval
        for m in d.get("measurements", []):
            set_measurement(p, user, m)
        for state_name, values in d.get("state_edits", {}).items():
            state = p.states.get(name=state_name)
            set_fields(state, values, STATE_EDITABLE)
            state.save()
        set_fields(p, d.get("protein_edits", {}), PROTEIN_EDITABLE)
        if doi := d.get("primary_doi"):
            # same call the public form makes (looks the DOI up on first use)
            old, new = p.primary_reference, Reference.objects.get_or_create(doi=doi.lower())[0]
            if old != new:
                p.primary_reference = new
                if old:
                    p.references.add(old)
                Lineage.objects.filter(protein=p, reference=old).update(reference=new)
        if mutation := d.get("lineage_mutation"):
            lineage = p.lineage
            if str(lineage.parent.protein.seq.mutate(mutation)) != str(p.seq):
                raise Skip(f"{lineage.parent.protein} + {mutation} does not give this sequence")
            lineage.mutation = mutation
            lineage.save()
        for doi in d.get("remove_references", []):
            ref = p.references.get(doi=doi.lower())
            if ref.id == p.primary_reference_id:
                raise Skip(f"{doi} is the primary reference")
            p.references.remove(ref)
        if approve:
            p.status = "approved"
        p.save()


def reject_protein(p, user, reason):
    approved = p.last_approved_version()
    if approved is None:
        # only hide genuinely new submissions: `status_changed` still equal to `created`
        # means it never left "pending" (see never_approved in fetch_pending.py)
        if (p.status_changed - p.created) > timedelta(minutes=1):
            raise Skip("established protein with no approved snapshot to revert to")
        with reversion.create_revision():
            reversion.set_user(user)
            reversion.set_comment(f"{user} rejected submission (hidden) {TAG} {reason}")
            p.status = "hidden"
            p.save()
        return "hidden"

    rev_ids = p.versions.values_list("revision_id", flat=True)
    later = Revision.objects.filter(id__in=rev_ids, id__gt=approved.revision_id)
    if staff := later.filter(user__is_staff=True).first():
        raise Skip(f"revision {staff.id} by staff user {staff.user} would be reverted too")

    revision = approved.revision
    revision.revert(delete=True)
    p.refresh_from_db()
    for state in p.states.all():
        state.save()  # revert writes raw rows; recompute derived fields (hex, brightness)
    if p.status != "approved":
        raise Skip(f"status is {p.status!r} after revert to revision {revision.id}")
    with reversion.create_revision():
        reversion.set_user(user)
        reversion.set_comment(f"Reverted to revision dated {revision.date_created} {TAG} {reason}")
        p.save()
    return f"reverted to revision {revision.id}"


def apply_one(d, user):
    reason = d.get("reason", "")
    if d["kind"] == "spectrum":
        sp = Spectrum.objects.all_objects().get(id=d["id"])
        if sp.status != "pending":
            raise Skip(f"status is {sp.status!r}, not pending")
        status = {"approve": Spectrum.STATUS.approved, "reject": Spectrum.STATUS.rejected}
        # save() rather than update(), so the post_save cache invalidation fires
        sp.status = status[d["action"]]
        sp.save()
        return {"result": sp.status, "uncache": getattr(sp.owner_fluor, "owner_slug", None)}

    p = Protein.objects.get(slug=d["slug"])
    if p.status != "pending" and d["action"] != "correct":
        raise Skip(f"status is {p.status!r}, not pending")
    if str(p.modified) != d["expect_modified"]:
        raise Skip(f"modified since fetch ({p.modified})")

    before = snapshot(p)
    approved = p.last_approved_version()
    allowed = pending_fields(p, approved) if d["action"] == "reject" and approved else None
    if d["action"] == "approve":
        approve_protein(p, user, reason, d)
        result = "approved"
    elif d["action"] == "correct":
        approve_protein(p, user, reason, d, approve=False)
        result = f"corrected ({p.status})"
    elif d["action"] == "reject":
        result = reject_protein(p, user, reason)
    else:
        raise Skip(f"unknown action {d['action']!r}")
    p = Protein.objects.get(id=p.id)
    after = snapshot(p)
    diff = {
        k: [before.get(k), after.get(k)]
        for k in sorted(before.keys() | after.keys())
        if before.get(k) != after.get(k)
    }
    if allowed and (bad := collateral(diff, allowed)):
        raise Skip(f"revert would also change {bad}: { ({k: diff[k] for k in bad}) }")
    return {"result": result, "data_diff": diff, "uncache": p.slug}


def uncache(slug):
    # best effort: only hits the anonymous-visitor cache key; pages expire in 30 min anyway
    request = HttpRequest()
    request.META = {
        "HTTP_HOST": "www.fpbase.org",
        "SERVER_PORT": "443",
        "HTTP_X_FORWARDED_PROTO": "https",
    }
    with contextlib.suppress(Exception):
        uncache_protein_page(slug, request)


commit = bool(PARAMS.get("commit"))
moderator = get_user_model().objects.get(username=PARAMS["moderator"], is_staff=True)
results = []

with transaction.atomic():
    for d in PARAMS["decisions"]:
        ident = {k: d[k] for k in ("kind", "slug", "id", "action") if k in d}
        try:
            # savepoint per decision: a failed revert must not poison the batch
            with transaction.atomic():
                results.append({**ident, "ok": True, **apply_one(d, moderator)})
        except Exception as e:
            results.append({**ident, "ok": False, "error": f"{type(e).__name__}: {e}"})
    if not commit:
        transaction.set_rollback(True)

if commit:
    for r in results:
        if r.get("uncache"):
            uncache(r["uncache"])

emit({"committed": commit, "moderator": moderator.username, "results": results})
