# ruff: noqa: F821  (PARAMS and emit are injected by remote.py / _bootstrap.py)
# READ-ONLY.  Dumps pending proteins/spectra as JSON.  Runs remotely (see remote.py).
# PARAMS: kind ("proteins"|"spectra"|"all"), limit, offset, slugs, summary
from datetime import timedelta

from proteins.models import Protein, Spectrum

# same keys that views.protein_history ignores
IGNORE_KEYS = [
    "modified",
    "created_by_id",
    "status_changed",
    "emhex",
    "exhex",
    "updated_by_id",
    "status",
    "seq_comment",
]


def fields(obj, skip=()):
    return {
        f.name: f.value_from_object(obj) for f in obj._meta.concrete_fields if f.name not in skip
    }


def user_info(user):
    if user is None:
        return None
    return {
        "username": user.username,
        "is_staff": user.is_staff,
        "date_joined": user.date_joined,
        "n_approved_proteins": Protein.objects.filter(created_by=user, status="approved").count(),
    }


def ref_info(ref):
    if ref is None:
        return None
    return {k: getattr(ref, k) for k in ("id", "doi", "pmid", "citation", "title", "year")}


def lineage_info(lin):
    return {
        "parent": lin.parent and lin.parent.protein.slug,
        "mutation": str(lin.mutation),
        "reference": ref_info(lin.reference),
    }


def never_approved(p):
    # `status_changed` only moves when status changes, so if it still equals `created` the
    # protein has never left "pending".  (Don't use reversion for this: admin bulk-approval
    # uses queryset.update(), so most approved proteins have no "approved" version.)
    unchanged = (p.status_changed - p.created) < timedelta(minutes=1)
    return unchanged and p.last_approved_version() is None


def protein_summary(p):
    return {
        "slug": p.slug,
        "name": p.name,
        "modified": p.modified,
        "is_new": never_approved(p),
        "updated_by": p.updated_by.username if p.updated_by else None,
        "primary_doi": p.primary_reference.doi if p.primary_reference else None,
    }


def protein_info(p):
    approved = p.last_approved_version()
    # no approved snapshot: fall back to revisions made since the record went pending
    since = p.created if never_approved(p) else p.status_changed - timedelta(minutes=1)
    pending_revisions = [
        {
            "id": rev.id,
            "date": rev.date_created,
            "user": rev.user and rev.user.username,
            "user_is_staff": bool(rev.user and rev.user.is_staff),
            "comment": rev.comment,
            "changes": dict(changes),
        }
        for rev, changes in p.history(IGNORE_KEYS).items()
        if (rev.id > approved.revision_id if approved else rev.date_created >= since)
    ]
    return {
        **protein_summary(p),
        "url": f"https://www.fpbase.org{p.get_absolute_url()}",
        "created_by": user_info(p.created_by),
        "last_approved_revision": approved
        and {
            "id": approved.revision_id,
            "date": approved.revision.date_created,
        },
        "n_versions": p.versions.count(),
        "pending_revisions": pending_revisions,
        "fields": fields(p),
        "parent_organism": str(p.parent_organism) if p.parent_organism else None,
        "primary_reference": ref_info(p.primary_reference),
        "references": [ref_info(r) for r in p.references.all()],
        "states": [fields(s) for s in p.states.all()],
        "transitions": [str(t) for t in p.transitions.all()],
        "lineage": lineage_info(p.lineage) if hasattr(p, "lineage") else None,
        "excerpts": [
            {
                "doi": e.reference.doi if e.reference else None,
                "content": e.content,
                "status": e.status,
                "created": e.created,
                "created_by": e.created_by and e.created_by.username,
            }
            for e in p.excerpts.all()
        ],
    }


def spectrum_info(sp):
    owner = sp.owner
    fluor = sp.owner_fluor
    data = sp.data or []
    return {
        **fields(sp, skip=("y_values",)),
        "owner": str(owner),
        "owner_type": type(owner).__name__,
        "owner_slug": getattr(fluor, "owner_slug", None),  # protein/dye slug
        "owner_ex_max": getattr(fluor, "ex_max", None),
        "owner_em_max": getattr(fluor, "em_max", None),
        "owner_approved_spectra": [
            f"{s.subtype} (peak {s.peak_wave})" for s in fluor.spectra.all()
        ]
        if fluor
        else None,
        "reference": ref_info(sp.reference),
        "created_by": user_info(sp.created_by),
        "admin_url": f"https://www.fpbase.org/admin/proteins/spectrum/{sp.id}/change/",
        # coarse trace, enough to sanity check shape
        "trace_10nm": [(x, round(y, 3)) for x, y in data if x % 10 == 0],
    }


def safely(func, obj, ident):
    # old reversion data can fail to deserialize; don't let one record kill the batch
    try:
        return func(obj)
    except Exception as e:
        return {**ident, "error": f"{type(e).__name__}: {e}"}


kind = PARAMS.get("kind", "all")
start = PARAMS.get("offset", 0)
stop = start + PARAMS["limit"] if PARAMS.get("limit") else None
out = {}

if kind in ("proteins", "all"):
    qs = (
        Protein.objects.filter(status="pending")
        .select_related("primary_reference", "updated_by", "created_by", "parent_organism")
        .order_by("-modified")
    )
    out["n_pending_proteins"] = qs.count()
    if PARAMS.get("slugs"):
        qs = qs.filter(slug__in=PARAMS["slugs"])
    func = protein_summary if PARAMS.get("summary") else protein_info
    out["proteins"] = [safely(func, p, {"slug": p.slug}) for p in qs[start:stop]]

if kind in ("spectra", "all") and not PARAMS.get("slugs"):
    qs = (
        Spectrum.objects.all_objects()
        .filter(status="pending")
        .select_related(
            "owner_fluor", "owner_filter", "owner_light", "owner_camera", "reference", "created_by"
        )
        .order_by("-created")
    )
    out["n_pending_spectra"] = qs.count()
    if not PARAMS.get("summary"):
        out["spectra"] = [safely(spectrum_info, s, {"id": s.id}) for s in qs[start:stop]]

emit(out)
