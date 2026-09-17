"""Payload for the client-side autocomplete search (replaces the Algolia indices)."""

from __future__ import annotations

import hashlib
import html
import json
import logging
import math
from collections import defaultdict
from typing import Any

from django.core.cache import cache
from django.db.models import Count, Prefetch
from django.urls import reverse
from django.utils.html import strip_tags

from favit.models import Favorite
from fpbase.cache_utils import get_model_version
from proteins.extrest.ga import cached_ga_popular, cached_ga_spectra_views
from proteins.models import Dye, DyeState, Organism, Protein, Spectrum, State
from proteins.util.helpers import get_color_group
from references.models import Reference

logger = logging.getLogger(__name__)

SEARCH_INDEX_VERSION = 2
CACHE_TTL = 60 * 60 * 24  # also bounds how stale popularity can get

# Popularity is log-scaled into [0, 1]: the most popular protein is 1, and each
# e-fold below it loses 1/POP_EFOLDS. ~11 e-folds spans the most-viewed protein
# (~60K views/yr) down to a single view.
POP_EFOLDS = 11.0
# Weight of a protein's share of all favorites, relative to its share of page views
FAVE_WEIGHT = 0.1


def _log_popularity(raw: dict[int, float]) -> dict[int, float]:
    top = max(raw.values(), default=0)
    if top <= 0:
        return {}
    return {
        k: round(max(0.0, 1 + math.log(v / top) / POP_EFOLDS), 3) for k, v in raw.items() if v > 0
    }


def _ga_view_shares() -> dict[str, float]:
    """Protein page view share (%) over the last year, keyed by slug."""
    try:
        return {slug: pct for slug, _name, pct in cached_ga_popular()["year"]}
    except Exception:
        # no credentials (dev/test) or GA outage: fall back to favorites only
        logger.warning("Could not fetch Google Analytics popularity", exc_info=True)
        return {}


def _ga_spectrum_views() -> dict[int, int]:
    """Spectra viewer page views over the last year, keyed by spectrum ID."""
    try:
        return cached_ga_spectra_views()
    except Exception:
        logger.warning("Could not fetch Google Analytics spectra views", exc_info=True)
        return {}


def protein_popularity(proteins: list[Protein]) -> dict[int, float]:
    views = _ga_view_shares()
    faves = dict(
        Favorite.objects.for_model(Protein)
        .order_by()  # Meta.ordering would otherwise leak into the GROUP BY
        .values("target_object_id")
        .annotate(n=Count("id"))
        .values_list("target_object_id", "n")
    )
    total_faves = sum(faves.values()) or 1
    raw = {
        p.id: views.get(p.slug, 0) + FAVE_WEIGHT * 100 * faves.get(p.id, 0) / total_faves
        for p in proteins
    }
    return _log_popularity(raw)


def _compact(d: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in d.items() if v not in (None, "", [], False)}


def _color(ex: float | None, em: float | None) -> str | None:
    """Same color group as `Protein.color`, without the default_state query."""
    if ex is None or em is None or not (group := get_color_group(ex, em)):
        return None
    return group[0]


def _protein_records(proteins: list[Protein], pop: dict[int, float]) -> list[dict]:
    states: dict[int, list[tuple]] = defaultdict(list)
    state_exem: dict[int, tuple] = {}
    for sid, pid, ex, em in State.objects.order_by("id").values_list(
        "id", "protein_id", "ex_max", "em_max"
    ):
        states[pid].append((ex, em))
        state_exem[sid] = (ex, em)
    with_spectra = set(
        Spectrum.objects.filter(owner_fluor__state__isnull=False, subtype__in=("ex", "em", "ab"))
        .values_list("owner_fluor__state__protein_id", flat=True)
        .distinct()
    )
    records = []
    for p in proteins:
        ex = [s[0] for s in states[p.id]]
        em = [s[1] for s in states[p.id]]
        records.append(
            _compact(
                {
                    "name": html.unescape(p.name),
                    "slug": p.slug,
                    "url": p.get_absolute_url(),
                    "aliases": [html.unescape(a) for a in p.aliases or []],
                    "uuid": p.uuid,
                    "pdb": p.pdb or [],
                    "genbank": p.genbank,
                    "uniprot": p.uniprot,
                    "ipg_id": p.ipg_id,
                    "color": _color(*state_exem.get(p.default_state_id, (None, None))),
                    "switch": p.switchType(),
                    "ex": ex[0] if len(ex) == 1 else ex,
                    "em": em[0] if len(em) == 1 else em,
                    "spectra": p.id in with_spectra,
                    "p": pop.get(p.id, 0),
                }
            )
        )
    return records


def _reference_records(pop: dict[int, float]) -> list[dict]:
    visible = Protein.objects.exclude(status=Protein.STATUS.hidden).only("id", "name")
    refs = Reference.objects.only(
        "id", "doi", "pmid", "title", "citation", "year"
    ).prefetch_related(
        Prefetch(
            "primary_proteins",
            queryset=visible.only("id", "name", "primary_reference_id"),
            to_attr="_primary",
        ),
        Prefetch("proteins", queryset=visible, to_attr="_all"),
    )
    records = []
    for ref in refs:
        primary_ids = {p.id for p in ref._primary}
        secondary = [p for p in ref._all if p.id not in primary_ids]
        # a paper is as popular as the proteins it introduced (discounted for mentions)
        ref_pop = max(
            [pop.get(p.id, 0) for p in ref._primary]
            + [0.5 * pop.get(p.id, 0) for p in secondary]
            + [0]
        )
        records.append(
            _compact(
                {
                    "citation": ref.citation,
                    "title": " ".join(strip_tags(ref.title).split()),
                    "doi": ref.doi,
                    "pmid": ref.pmid,
                    "year": ref.year,
                    "url": ref.get_absolute_url(),
                    "primary": [html.unescape(p.name) for p in ref._primary],
                    "secondary": [html.unescape(p.name) for p in secondary],
                    "p": round(ref_pop, 3),
                }
            )
        )
    return records


def _organism_records() -> list[dict]:
    orgs = list(
        Organism.objects.annotate(n=Count("proteins")).values_list("id", "scientific_name", "n")
    )
    pop = _log_popularity({oid: n for oid, _, n in orgs})
    return [
        _compact(
            {
                "name": name,
                "url": reverse("proteins:organism-detail", args=[oid]),
                "p": pop.get(oid, 0),
            }
        )
        for oid, name, n in orgs
    ]


def _dye_records() -> list[dict]:
    """Dyes have no detail page, so they link to their spectra in the spectra viewer."""
    spectra: dict[int, list[int]] = defaultdict(list)
    for spectrum_id, state_id in (
        Spectrum.objects.filter(
            owner_fluor__dyestate__isnull=False, subtype__in=("ex", "em", "ab")
        )
        .order_by("id")
        .values_list("id", "owner_fluor_id")
    ):
        spectra[state_id].append(spectrum_id)
    # use the default state if it has spectra, otherwise the first state that does
    default_state = dict(Dye.objects.values_list("id", "default_state_id"))
    chosen: dict[int, tuple] = {}
    for row in DyeState.objects.order_by("id").values_list(
        "id", "dye_id", "dye__name", "ex_max", "em_max", "emhex"
    ):
        state_id, dye_id = row[:2]
        if spectra[state_id] and (dye_id not in chosen or state_id == default_state[dye_id]):
            chosen[dye_id] = row
    views = _ga_spectrum_views()
    pop = _log_popularity(
        {dye_id: max(views.get(s, 0) for s in spectra[row[0]]) for dye_id, row in chosen.items()}
    )
    viewer = reverse("proteins:spectra")
    return [
        _compact(
            {
                "name": html.unescape(name),
                "url": f"{viewer}?s={','.join(map(str, spectra[state_id]))}",
                "ex": ex,
                "em": em,
                "color": emhex if em else None,  # emhex is a placeholder without em_max
                "p": pop.get(dye_id, 0),
            }
        )
        for dye_id, (state_id, _, name, ex, em, emhex) in chosen.items()
    ]


def build_search_index() -> dict[str, Any]:
    proteins = list(
        Protein.objects.exclude(status=Protein.STATUS.hidden).only(
            "id", "name", "slug", "aliases", "uuid", "pdb", "genbank", "uniprot", "ipg_id",
            "default_state_id", "switch_type",
        )
    )  # fmt: skip
    pop = protein_popularity(proteins)
    return {
        "version": SEARCH_INDEX_VERSION,
        "proteins": _protein_records(proteins, pop),
        "references": _reference_records(pop),
        "organisms": _organism_records(),
        "dyes": _dye_records(),
    }


def get_search_index() -> tuple[bytes, str]:
    """Return the serialized search index and its ETag, rebuilding if models changed."""
    version = get_model_version(Protein, State, Reference, Organism, Dye, DyeState, Spectrum)
    key = f"search_index:{SEARCH_INDEX_VERSION}:{version}"
    if (hit := cache.get(key)) is None:
        data = json.dumps(build_search_index(), separators=(",", ":")).encode()
        hit = (data, f'W/"{hashlib.blake2b(data, digest_size=16).hexdigest()}"')
        cache.set(key, hit, CACHE_TTL)
    return hit
