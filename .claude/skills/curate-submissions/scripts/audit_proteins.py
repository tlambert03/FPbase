# ruff: noqa: F821  (PARAMS and emit are injected by remote.py / _bootstrap.py)
# READ-ONLY.  Machine-checkable health of given proteins (any status): what is missing,
# and does what is there agree with itself?  No literature, no judgement.
from proteins.models import Protein, Spectrum

CORE_STATE_FIELDS = ("ex_max", "em_max", "ext_coeff", "qy", "pka", "lifetime", "maturation")


def lineage_matches_seq(p):
    lin = getattr(p, "lineage", None)
    if not (lin and lin.parent and lin.mutation and p.seq and lin.parent.protein.seq):
        return None
    try:
        return str(lin.parent.protein.seq.mutate(lin.mutation)) == str(p.seq)
    except Exception as e:
        return f"error: {type(e).__name__}"


def audit(p):
    state = p.default_state
    spectra = Spectrum.objects.filter(owner_fluor__in=p.states.values("fluorstate_ptr"))
    return {
        "slug": p.slug,
        "name": p.name,
        "status": p.status,
        "seq": str(p.seq) if p.seq else None,
        "seq_validated": p.seq_validated,
        "genbank": p.genbank,
        "uniprot": p.uniprot,
        "pdb": p.pdb,
        "lineage_matches_seq": lineage_matches_seq(p),
        "has_lineage": hasattr(p, "lineage"),
        "primary_doi": p.primary_reference and p.primary_reference.doi,
        "agg": p.agg,
        "organism": p.parent_organism_id,
        "n_states": p.states.count(),
        "missing": [f for f in CORE_STATE_FIELDS if state is None or getattr(state, f) is None],
        "is_dark": bool(state and state.is_dark),
        "spectra": sorted(spectra.values_list("subtype", flat=True)),
    }


qs = Protein.objects.filter(slug__in=PARAMS["slugs"]).select_related(
    "default_state", "primary_reference", "lineage__parent__protein"
)
emit([audit(p) for p in qs])
