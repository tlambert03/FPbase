"""Tests for collection views to prevent N+1 query regressions."""

from __future__ import annotations

from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory
from proteins.models import ProteinCollection


class CollectionDetailViewTests(TestCase):
    """Test the CollectionDetail view for performance."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests."""
        cls.collection = ProteinCollection.objects.create(
            name="Test Collection",
            description="A test collection",
            private=False,
        )

        for _ in range(10):
            protein = ProteinFactory()
            StateFactory(protein=protein, name="state1")
            cls.collection.proteins.add(protein)

    @override_settings(DEBUG=True)
    def test_collection_json_export_query_count(self):
        """
        Test that collection JSON export doesn't generate excessive queries.

        This prevents N+1 query regressions. With proper prefetch_related usage,
        the query count should be constant (< 10) regardless of collection size.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.get(f"/collection/{self.collection.pk}/?format=json")

        self.assertEqual(response.status_code, 200)

        query_count = len(context.captured_queries)
        self.assertLessEqual(
            query_count,
            8,
            f"Collection JSON export generated {query_count} queries. "
            f"Expected ≤8 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    @override_settings(DEBUG=True)
    def test_collection_csv_export_query_count(self):
        """Test that collection CSV export doesn't generate excessive queries."""
        with CaptureQueriesContext(connection) as context:
            response = self.client.get(f"/collection/{self.collection.pk}/?format=csv")

        self.assertEqual(response.status_code, 200)

        query_count = len(context.captured_queries)
        self.assertLessEqual(
            query_count,
            8,
            f"Collection CSV export generated {query_count} queries. "
            f"Expected ≤8 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )
