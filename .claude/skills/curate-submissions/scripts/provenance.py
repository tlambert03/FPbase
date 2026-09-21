# ruff: noqa: F821  (PARAMS and emit are injected by remote.py / _bootstrap.py)
# READ-ONLY.  Where did a record's values come from?  For one protein: every change to
# each field (who, when, old -> new), the spectra attached to each state (raw peak, smoothed
# peak, scale factor, who uploaded them), and who the editors are.
# PARAMS: slug, around (optional list of wavelengths to print spectrum data around)
import json

from reversion.models import Version

from proteins.models import Protein, Spectrum, State

FLUOR_STATE = State._meta.get_parent_list()[0]
PROTEIN_FIELDS = (
    "name", "aliases", "seq", "pdb", "genbank", "uniprot", "ipg_id", "agg", "switch_type",
    "cofactor", "chromophore", "parent_organism", "primary_reference", "status",
)  # fmt: skip
STATE_FIELDS = (
    "name", "ex_max", "em_max", "ext_coeff", "qy", "pka", "lifetime", "maturation",
    "twop_ex_max", "twop_peak_gm", "twop_qy", "is_dark",
)  # fmt: skip


def short(v):
    s = str(v)
    return s if len(s) <= 40 else f"{s[:18]}…{s[-12:]} ({len(s)} chars)"


def field_history(versions_by_model, fields, editors):
    """Changes to `fields` across snapshots: [{date, user, comment, changes: {f: [old, new]}}]."""
    snaps = {}  # revision id -> (revision, merged fields)
    for model, pk in versions_by_model:
        qs = Version.objects.get_for_object_reference(model, pk).select_related("revision__user")
        for v in qs:
            rev, merged = snaps.setdefault(v.revision_id, (v.revision, {}))
            merged.update(json.loads(v.serialized_data)[0]["fields"])
    out, prev = [], None
    for _, (rev, now) in sorted(snaps.items()):
        if rev.user:
            u = rev.user
            editors[u.username] = {
                "email": u.email,
                "name": u.get_full_name(),
                "joined": u.date_joined,
                "is_staff": u.is_staff,
            }
        changes = {
            f: [short(prev.get(f)) if prev else None, short(now[f])]
            for f in fields
            if f in now
            and (prev is None or prev.get(f) != now[f])
            and (prev or now[f] not in (None, "", []))
        }
        if changes:
            out.append(
                {
                    "date": rev.date_created,
                    "user": rev.user and rev.user.username,
                    "first_snapshot": prev is None,
                    "changes": changes,
                    "comment": (rev.comment or "")[:200],
                }
            )
        prev = now
    return out


def measurement_info(m):
    keep = ("ex_max", "em_max", "ext_coeff", "qy", "pka", "lifetime", "twop_ex_max",
            "twop_peak_gm", "twop_qy", "is_dark")  # fmt: skip
    return {
        "id": m.id,
        "doi": m.reference and m.reference.doi,
        "values": {k: getattr(m, k) for k in keep if getattr(m, k) not in (None, False)},
        "conditions": m.conditions,
        "created_by": m.created_by and m.created_by.username,
        "modified": m.modified,
    }


def smoothed_peak(data, near, window=11, reach=40):
    """Maximum of a moving average within `reach` nm of the raw peak `near`.

    Robust to one-point spikes.  Local on purpose: a 2P spectrum's global smooth maximum
    is usually its short-wavelength edge, which is not the peak anyone means.
    """
    xs, ys = [x for x, _ in data], [y for _, y in data]
    if near is None or len(ys) < window:
        return None
    half = window // 2
    means = {
        xs[i]: sum(ys[i - half : i + half + 1]) / window
        for i in range(half, len(ys) - half)
        if abs(xs[i] - near) <= reach
    }
    return max(means, key=means.get) if means else None


def plateau(data, near, level=0.98):
    """Contiguous wavelength range around the raw peak where y >= level * y(peak).

    A wide plateau means the data cannot say where "the" peak is; any value inside it is
    consistent with the spectrum.
    """
    lookup = dict(data)
    if near not in lookup:
        return None
    floor = level * lookup[near]
    lo = hi = near
    while lookup.get(lo - 1, 0) >= floor:
        lo -= 1
    while lookup.get(hi + 1, 0) >= floor:
        hi += 1
    return [lo, hi]


def spectrum_info(sp, around):
    data = sp.data or []
    lookup = dict(data)
    return {
        "id": sp.id,
        "subtype": sp.subtype,
        "status": sp.status,
        "peak_wave": sp.peak_wave,  # raw argmax, what the site displays
        "smoothed_peak": smoothed_peak(data, sp.peak_wave),
        "plateau_98pct": plateau(data, sp.peak_wave),
        "scale_factor": sp.scale_factor,
        "range": [sp.min_wave, sp.max_wave],
        "created": sp.created,
        "created_by": sp.created_by and sp.created_by.username,
        "reference": sp.reference and sp.reference.doi,
        "source": sp.source,
        "windows": {
            w: {x: round(lookup[x], 4) for x in range(w - 16, w + 17, 2) if x in lookup}
            for w in around
            if sp.min_wave <= w <= sp.max_wave
        },
    }


p = Protein.objects.get(slug=PARAMS["slug"])
around = PARAMS.get("around") or []
editors = {}
out = {
    "slug": p.slug,
    "name": p.name,
    "status": p.status,
    "modified": p.modified,  # pass as expect_modified to `apply`
    "created": p.created,
    "status_changed": p.status_changed,
    "protein_history": field_history([(Protein, p.pk)], PROTEIN_FIELDS, editors),
    "states": [
        {
            "name": s.name,
            "now": {f: getattr(s, f) for f in STATE_FIELDS},
            # spectral fields live on State in old snapshots and on FluorState in newer ones
            "history": field_history([(State, s.pk), (FLUOR_STATE, s.pk)], STATE_FIELDS, editors),
            "measurements": [
                measurement_info(m)
                for m in s.measurements.select_related("reference", "created_by")
            ],
            "spectra": [
                spectrum_info(sp, around)
                for sp in Spectrum.objects.all_objects().filter(owner_fluor_id=s.pk)
            ],
        }
        for s in p.states.all()
    ],
}
out["editors"] = editors
emit(out)
