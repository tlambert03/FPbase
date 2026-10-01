"""What the public may not see: proteins with status "hidden", and what belongs to them."""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

from proteins.models import OSERMeasurement, Protein, Spectrum, State

if TYPE_CHECKING:
    from django.http import HttpRequest


class HiddenIds(NamedTuple):
    proteins: frozenset[int]
    states: frozenset[int]
    spectra: frozenset[int]
    oser_measurements: frozenset[int]


def hidden_ids(request: HttpRequest) -> HiddenIds:
    """Ids of the hidden proteins and of their states, spectra and OSER measurements.

    For lists that were already fetched (prefetched relations in GraphQL): filtering
    their querysets instead would discard the prefetch and run a query per row.
    Looked up once per request.
    """
    if (ids := getattr(request, "_hidden_ids", None)) is None:
        hidden = Protein.objects.filter(status=Protein.STATUS.hidden)
        proteins = frozenset(hidden.values_list("id", flat=True))
        states = State.objects.filter(protein__in=proteins).values_list("id", flat=True)
        spectra = Spectrum.objects.all_objects().filter(owner_fluor__in=states)
        oser = OSERMeasurement.objects.filter(protein__in=proteins)
        ids = HiddenIds(
            proteins,
            frozenset(states),
            frozenset(spectra.values_list("id", flat=True)),
            frozenset(oser.values_list("id", flat=True)),
        )
        request._hidden_ids = ids  # pyright: ignore[reportAttributeAccessIssue]
    return ids
