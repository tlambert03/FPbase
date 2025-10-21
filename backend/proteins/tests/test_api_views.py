"""
Tests for protein API views.

These tests ensure the API endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory


class ProteinListAPIViewTests(TestCase):
    """Test the ProteinListAPIView endpoint."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests in this class."""
        for _ in range(15):
            protein = ProteinFactory()
            StateFactory(protein=protein, name="state1")

    @override_settings(DEBUG=True)
    def test_protein_list_api_query_count(self):
        """
        Test that the protein list API doesn't generate excessive queries.

        This prevents N+1 query regressions. With proper prefetch_related usage,
        the query count should be constant (4-5 queries) regardless of protein count.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.get("/api/proteins/", headers={"accept": "application/json"})

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data), 15)

        query_count = len(context.captured_queries)
        self.assertLess(
            query_count,
            5,
            f"API endpoint generated {query_count} queries. "
            f"Expected 4 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    def test_protein_list_api_returns_expected_fields(self):
        """Test that the API returns the expected fields for each protein."""
        response = self.client.get("/api/proteins/", headers={"accept": "application/json"})

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data), 15)

        if data:
            protein = data[0]
            expected_fields = {"uuid", "name", "slug", "states", "transitions"}
            self.assertTrue(
                expected_fields.issubset(protein.keys()),
                f"Expected fields {expected_fields} not all present in {protein.keys()}",
            )

            self.assertGreater(len(protein.get("states", [])), 0, "States should be included")
