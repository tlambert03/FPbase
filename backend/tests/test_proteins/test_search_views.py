"""Tests for protein search views to prevent N+1 query regressions."""

from __future__ import annotations

from django.db import connection
from django.db.models import Count
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory
from proteins.filters import ProteinFilter
from proteins.models import Protein


class ProteinSearchViewTests(TestCase):
    """Test the protein_search queryset optimization for performance."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests."""
        for i in range(10):
            protein = ProteinFactory(name=f"TestProtein{i}")
            StateFactory(protein=protein, name="state1")
            StateFactory(protein=protein, name="state2")

    @override_settings(DEBUG=True)
    def test_protein_search_queryset_optimization(self):
        """
        Test that protein search queryset doesn't generate excessive queries.

        This prevents N+1 query regressions. With proper prefetch_related usage,
        accessing protein attributes should not trigger additional queries.
        """
        f = ProteinFilter(
            {"name__icontains": "TestProtein"},
            queryset=Protein.visible.annotate(nstates=Count("states"))
            .select_related("default_state", "primary_reference")
            .prefetch_related("states__spectra", "transitions")
            .order_by("default_state__em_max"),
        )

        with CaptureQueriesContext(connection) as context:
            proteins = list(f.qs)
            for protein in proteins:
                _ = protein.name
                _ = protein.default_state
                _ = protein.slug
                if protein.default_state:
                    _ = list(protein.states.all())
                    for state in protein.states.all():
                        _ = list(state.spectra.all())

        query_count = len(context.captured_queries)
        self.assertLessEqual(
            query_count,
            10,
            f"Protein search queryset generated {query_count} queries. "
            f"Expected <= 10 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )
