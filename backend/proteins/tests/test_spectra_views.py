"""Tests for spectra views to prevent N+1 query regressions."""

from __future__ import annotations

from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import SpectrumFactory, StateFactory
from proteins.models import Spectrum


class SpectraCSVViewTests(TestCase):
    """Test the spectra_csv view for performance."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests."""
        cls.spectra_list = []
        for i in range(10):
            state = StateFactory()
            spectrum = SpectrumFactory(
                owner_state=state,
                subtype=Spectrum.EX if i % 2 == 0 else Spectrum.EM,
                category=Spectrum.PROTEIN,
            )
            cls.spectra_list.append(spectrum)

    @override_settings(DEBUG=True)
    def test_spectra_csv_query_count(self):
        """
        Test that spectra CSV export doesn't generate excessive queries.

        This prevents N+1 query regressions. With proper select_related usage,
        the query count should be constant (< 5) regardless of spectrum count.
        """
        spectrum_ids = ",".join(str(s.id) for s in self.spectra_list)

        with CaptureQueriesContext(connection) as context:
            response = self.client.get(f"/spectra_csv/?q={spectrum_ids}")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "text/csv")

        query_count = len(context.captured_queries)

        print(f"\n=== Query count: {query_count} ===")
        for i, query in enumerate(context.captured_queries, 1):
            print(f"{i}. {query['sql'][:150]}...")

        self.assertLessEqual(
            query_count,
            5,
            f"Spectra CSV export generated {query_count} queries. "
            f"Expected <= 5 queries with proper select_related. "
            f"This may indicate an N+1 query regression.",
        )
