"""Tests for the scripts in .claude/skills/curate-submissions (piped into `heroku run`)."""

from __future__ import annotations

import sys
from datetime import timedelta
from pathlib import Path

import pytest
import reversion
from django.contrib.auth import get_user_model
from reversion.models import Revision, Version

from proteins.factories import SpectrumFactory
from proteins.models import Lineage, Protein, Spectrum, State, StateTransition
from references.models import Reference

SCRIPTS = Path(__file__).parents[3] / ".claude/skills/curate-submissions/scripts"
User = get_user_model()


def run_script(name: str, params: dict) -> dict:
    out: list[dict] = []
    source = (SCRIPTS / name).read_text()
    exec(compile(source, name, "exec"), {"PARAMS": params, "emit": out.append})
    return out[0]


def apply(decision: dict, commit: bool = True) -> dict:
    params = {"decisions": [decision], "commit": commit, "moderator": "staff"}
    return run_script("apply_decisions.py", params)["results"][0]


@pytest.fixture
def staff(db) -> User:
    return User.objects.create_user(username="staff", is_staff=True)


@pytest.fixture
def submitter(db) -> User:
    return User.objects.create_user(username="submitter")


@pytest.fixture
def edited_protein(staff: User, submitter: User) -> Protein:
    """Approved protein with a pending edit (agg + ex_max) by a non-staff user."""
    with reversion.create_revision():
        reversion.set_user(staff)
        p = Protein.objects.create(name="EditedFP", agg="m", status="approved")
        State.objects.create(protein=p, name="default", ex_max=488, em_max=510)
    with reversion.create_revision():
        reversion.set_user(submitter)
        state = p.states.get()
        state.ex_max = 600
        state.save()
        p.agg = "d"
        p.status = "pending"
        p.save()
    return p


@pytest.fixture
def new_protein(staff: User, submitter: User) -> Protein:
    with reversion.create_revision():
        reversion.set_user(submitter)
        return Protein.objects.create(name="BrandNewFP", status="pending")


def decision(p: Protein, action: str) -> dict:
    p.refresh_from_db()
    return {
        "kind": "protein",
        "slug": p.slug,
        "action": action,
        "expect_modified": str(p.modified),
    }


def test_fetch_pending(edited_protein: Protein, new_protein: Protein) -> None:
    out = run_script("fetch_pending.py", {"kind": "proteins"})
    assert out["n_pending_proteins"] == 2
    by_slug = {p["slug"]: p for p in out["proteins"]}
    assert not any("error" in p for p in by_slug.values())
    assert by_slug[new_protein.slug]["is_new"]
    edited = by_slug[edited_protein.slug]
    assert not edited["is_new"]
    assert [r["user"] for r in edited["pending_revisions"]] == ["submitter"]
    assert edited["states"][0]["ex_max"] == 600


def test_approve_keeps_edit(edited_protein: Protein) -> None:
    result = apply(decision(edited_protein, "approve"))
    assert result["ok"], result
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.agg) == ("approved", "d")
    assert edited_protein.states.get().ex_max == 600


def test_approve_with_state_edits(edited_protein: Protein) -> None:
    edits = {"default": {"ext_coeff": 86100, "is_dark": True}}
    result = apply({**decision(edited_protein, "approve"), "state_edits": edits})
    assert result["ok"], result
    assert result["data_diff"]["state[default].ext_coeff"] == [None, 86100]
    state = edited_protein.states.get()
    assert (state.ext_coeff, state.is_dark, state.ex_max) == (86100, True, 600)

    bad = {"default": {"ext_coef": 1}}
    edited_protein.status = "pending"
    edited_protein.save()
    result = apply({**decision(edited_protein, "approve"), "state_edits": bad})
    assert not result["ok"]
    edited_protein.refresh_from_db()
    assert edited_protein.status == "pending"


def test_approve_undoing_the_edit(edited_protein: Protein) -> None:
    """ "undo edit": put the pre-pending values back and approve, without a reversion revert."""
    ref = Reference(doi="10.1234/wrong.paper", year=2020)
    ref.save(skipdoi=True)
    edited_protein.references.add(ref)
    result = apply(
        {
            **decision(edited_protein, "approve"),
            "protein_edits": {"agg": "m"},
            "state_edits": {"default": {"ex_max": 488}},
            "remove_references": [ref.doi],
        }
    )
    assert result["ok"], result
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.agg) == ("approved", "m")
    assert edited_protein.states.get().ex_max == 488
    assert not edited_protein.references.exists()

    edited_protein.status = "pending"
    edited_protein.save()
    result = apply({**decision(edited_protein, "approve"), "protein_edits": {"status": "hidden"}})
    assert not result["ok"] and "not an editable field" in result["error"]


def test_reject_edit_reverts(edited_protein: Protein) -> None:
    result = apply(decision(edited_protein, "reject"))
    assert result["ok"], result
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.agg) == ("approved", "m")
    assert edited_protein.states.get().ex_max == 488
    assert result["data_diff"]["protein.agg"] == ["d", "m"]


def test_reject_new_hides(new_protein: Protein) -> None:
    result = apply(decision(new_protein, "reject"))
    assert result["ok"], result
    new_protein.refresh_from_db()
    assert new_protein.status == "hidden"


def test_reject_never_hides_established_protein(staff: User, submitter: User) -> None:
    # approved without any reversion snapshot (like the admin bulk action), then edited
    p = Protein.objects.create(name="OldFP", status="approved")
    Protein.objects.filter(id=p.id).update(created=p.created - timedelta(days=900))
    p.refresh_from_db()
    with reversion.create_revision():
        reversion.set_user(submitter)
        p.status = "pending"
        p.save()
    out = run_script("fetch_pending.py", {"kind": "proteins"})
    assert not out["proteins"][0]["is_new"]
    result = apply(decision(p, "reject"))
    assert not result["ok"]
    p.refresh_from_db()
    assert p.status == "pending"


def test_dry_run_changes_nothing(edited_protein: Protein) -> None:
    result = apply(decision(edited_protein, "reject"), commit=False)
    assert result["ok"], result
    assert result["data_diff"]["protein.agg"] == ["d", "m"]
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.agg) == ("pending", "d")


def test_skips_if_modified_since_fetch(edited_protein: Protein) -> None:
    d = decision(edited_protein, "approve")
    edited_protein.save()  # bumps `modified`
    result = apply(d)
    assert not result["ok"]
    edited_protein.refresh_from_db()
    assert edited_protein.status == "pending"


def test_reject_refuses_to_revert_staff_edits(edited_protein: Protein, staff: User) -> None:
    with reversion.create_revision():
        reversion.set_user(staff)
        edited_protein.blurb = "staff fix"
        edited_protein.save()
    result = apply(decision(edited_protein, "reject"))
    assert not result["ok"]
    assert "staff" in result["error"]
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.blurb) == ("pending", "staff fix")


def test_reject_refuses_collateral_changes(edited_protein: Protein) -> None:
    # a change no revision knows about (stands in for schema drift since the approved revision)
    Protein.objects.filter(id=edited_protein.id).update(blurb="unversioned")
    result = apply(decision(edited_protein, "reject"))
    assert not result["ok"]
    assert "protein.blurb" in result["error"]
    edited_protein.refresh_from_db()
    assert (edited_protein.status, edited_protein.agg) == ("pending", "d")


def test_spectrum_fetch_and_approve(edited_protein: Protein) -> None:
    sp = SpectrumFactory(
        owner_fluor=edited_protein.states.get(),
        category=Spectrum.PROTEIN,
        subtype=Spectrum.EX,
        status=Spectrum.STATUS.pending,
    )
    out = run_script("fetch_pending.py", {"kind": "spectra"})
    assert out["n_pending_spectra"] == 1
    assert "error" not in out["spectra"][0], out["spectra"][0]
    assert out["spectra"][0]["owner_slug"] == edited_protein.slug

    result = apply({"kind": "spectrum", "id": sp.id, "action": "approve"})
    assert result["ok"], result
    assert Spectrum.objects.all_objects().get(id=sp.id).status == Spectrum.STATUS.approved


@pytest.fixture
def old_edited_protein(edited_protein: Protein) -> Protein:
    """`edited_protein`, but approved long ago (realistic: triage keys off timestamps)."""
    long_ago = edited_protein.created - timedelta(days=900)
    Protein.objects.filter(id=edited_protein.id).update(created=long_ago)
    State.objects.filter(protein=edited_protein).update(created=long_ago)
    first = Version.objects.get_for_object(edited_protein).last().revision
    Revision.objects.filter(id=first.id).update(date_created=long_ago)
    edited_protein.refresh_from_db()
    return edited_protein


def test_triage_reports_net_change(old_edited_protein: Protein, new_protein: Protein) -> None:
    rows = {r["slug"]: r for r in run_script("triage_pending.py", {})}
    assert rows[new_protein.slug]["kind"] == "new"
    edit = rows[old_edited_protein.slug]
    assert edit["kind"] == "edit" and edit["baseline"], edit
    assert edit["changes"]["protein.agg"] == ["m", "d"]
    assert edit["changes"]["state[default].ex_max"] == [488, 600]


def test_triage_no_net_change(old_edited_protein: Protein) -> None:
    # the submitter's values are put back (like the #422 restore): nothing left to review
    old_edited_protein.agg = "m"
    old_edited_protein.save()
    state = old_edited_protein.states.get()
    state.ex_max = 488
    state.save()
    (row,) = run_script("triage_pending.py", {})
    assert row.get("changes") == {}, row


def test_triage_sees_changed_and_deleted_transitions(
    old_edited_protein: Protein, staff: User
) -> None:
    p = old_edited_protein
    long_ago = p.created
    with reversion.create_revision():
        dark = State.objects.create(protein=p, name="dark", is_dark=True)
        tr = StateTransition.objects.create(
            protein=p, from_state=p.states.get(name="default"), to_state=dark, trans_wave=405
        )
    Revision.objects.filter(id=Revision.objects.latest("id").id).update(date_created=long_ago)
    State.objects.filter(id=dark.id).update(created=long_ago)
    StateTransition.objects.filter(id=tr.id).update(created=long_ago)

    with reversion.create_revision():
        tr.trans_wave = 488
        tr.save()
    (row,) = run_script("triage_pending.py", {})
    assert [v for k, v in row["changes"].items() if k.endswith(".trans_wave")] == [[405, 488]]

    tr.delete()
    (row,) = run_script("triage_pending.py", {})
    assert "transitions.deleted" in row["changes"]


def test_triage_deleted_ignores_other_proteins(old_edited_protein: Protein) -> None:
    """A state deleted from protein 61 must not show up as deleted from protein 6."""
    p = old_edited_protein
    with reversion.create_revision():
        other = Protein.objects.create(id=p.id * 10 + 1, name="OtherFP", status="approved")
        doomed = State.objects.create(protein=other, name="doomed", ex_max=400, em_max=450)
    Revision.objects.filter(id=Revision.objects.latest("id").id).update(date_created=p.created)
    doomed.delete()
    (row,) = run_script("triage_pending.py", {"slugs": [p.slug]})
    assert "states.deleted" not in row["changes"], row["changes"]


def test_triage_names_added_reference(old_edited_protein: Protein, submitter: User) -> None:
    ref = Reference(doi="10.1234/added.paper", year=2021)
    ref.save(skipdoi=True)
    with reversion.create_revision():
        reversion.set_user(submitter)
        old_edited_protein.references.add(ref)
        old_edited_protein.save()
    (row,) = run_script("triage_pending.py", {})
    (added,) = row["changes"]["references.added"][1]
    assert added["doi"] == ref.doi
    assert row["editors"] == ["submitter"]
    assert [c["doi"] for c in row["cited"]] == [ref.doi]


def test_lineage_fix_must_reproduce_sequence(staff: User, submitter: User) -> None:
    parent_seq = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT"
    parent = Protein.objects.create(name="ParentFP", seq=parent_seq, status="approved")
    Lineage.objects.create(protein=parent)
    child_seq = parent_seq.replace("MVSKGEE", "MVSRGDE")  # K4R/E6D
    with reversion.create_revision():
        reversion.set_user(submitter)
        child = Protein.objects.create(name="ChildFP", seq=child_seq, status="pending")
        Lineage.objects.create(protein=child, parent=parent.lineage, mutation="K4R")  # incomplete

    result = apply({**decision(child, "approve"), "lineage_mutation": "K4R/E6A"})  # wrong
    assert not result["ok"] and "does not give this sequence" in result["error"]
    child.refresh_from_db()
    assert child.status == "pending"

    result = apply(
        {
            **decision(child, "approve"),
            "lineage_mutation": "K4R/E6D",
            "protein_edits": {"name": "ChildFP2", "seq_validated": True},
        }
    )
    assert result["ok"], result
    child.refresh_from_db()
    assert (child.status, child.name, child.slug) == ("approved", "ChildFP2", "childfp2")
    assert child.seq_validated
    assert str(child.lineage.mutation) == "K4R/E6D"


def test_views_parse_and_rank() -> None:
    sys.path.insert(0, str(SCRIPTS))
    from views import annotate, parse_report  # local helper, not a Django module

    def row(path: str, n: int) -> dict:
        return {"dimension_values": [{"value": path}], "metric_values": [{"value": str(n)}]}

    report = {
        "rows": [
            row("/protein/mcherry/", 100),
            row("/protein/mCherry/", 5),  # same page, different case
            row("/protein/mcherry/history/", 50),  # not the detail page
            row("/protein/egfp/", 200),
        ]
    }
    views = parse_report(report)
    assert views == {"egfp": 200, "mcherry": 105}

    rows = annotate([{"slug": "mcherry"}, {"slug": "brand-new"}, {"slug": "egfp"}], views)
    assert [(r["slug"], r["views_365d"], r["rank"]) for r in rows] == [
        ("egfp", 200, 1),
        ("mcherry", 105, 2),
        ("brand-new", 0, None),
    ]


def test_audit_reports_gaps_and_consistency(old_edited_protein: Protein) -> None:
    (row,) = run_script("audit_proteins.py", {"slugs": [old_edited_protein.slug]})
    assert row["status"] == "pending"
    assert row["seq"] is None and not row["seq_validated"]
    assert "ext_coeff" in row["missing"] and "ex_max" not in row["missing"]
    assert row["lineage_matches_seq"] is None  # no lineage: nothing to check
