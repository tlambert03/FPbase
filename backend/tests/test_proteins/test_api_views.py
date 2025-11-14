"""
Tests for protein API views.

These tests ensure the API endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import ProteinFactory, StateFactory
from proteins.models import Protein


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
        self.assertLessEqual(
            query_count,
            6,
            f"API endpoint generated {query_count} queries. "
            f"Expected ≤6 queries with proper prefetch_related. "
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


class SpectraListAPIViewTests(TestCase):
    """Test the spectra-list endpoint with ETag caching."""

    @classmethod
    def setUpTestData(cls):
        """Create test proteins with states (which have spectra)."""
        cls.proteins = []
        for i in range(3):
            protein = ProteinFactory(name=f"TestProtein{i}")
            StateFactory(protein=protein, name="default")
            cls.proteins.append(protein)

    def test_spectra_list_etag_caching(self):
        """Test that ETag caching works correctly for spectra list.

        Verifies:
        1. Initial request returns 200 with ETag
        2. Second request with matching ETag returns 304
        3. After data changes, request with old ETag returns 200 with new data
        """
        # Initial request - should return 200 with ETag
        response1 = self.client.get("/api/spectra-list/")

        self.assertEqual(response1.status_code, 200)
        self.assertIn("ETag", response1)
        self.assertTrue(response1["ETag"].startswith('W/"'), "ETag should be a weak ETag")

        etag1 = response1["ETag"]
        data1 = response1.json()
        self.assertIn("data", data1)
        self.assertIn("spectra", data1["data"])
        names = {s["owner"]["name"] for s in data1["data"]["spectra"]}
        self.assertNotIn(
            "ModifiedProtein",
            names,
            "Should NOT find spectra with modified protein name",
        )

        # Second request with If-None-Match - should return 304
        response2 = self.client.get("/api/spectra-list/", headers={"If-None-Match": etag1})

        self.assertEqual(response2.status_code, 304)
        self.assertEqual(response2.content, b"")
        self.assertEqual(response2["ETag"], etag1, "304 response should include same ETag")

        # Modify a protein to invalidate cache
        protein = Protein.objects.first()
        protein.name = "ModifiedProtein"
        protein.save()

        # Third request with old ETag - should return 200 with new data
        response3 = self.client.get("/api/spectra-list/", headers={"If-None-Match": etag1})

        self.assertEqual(response3.status_code, 200, "Should return 200 after data changed")
        self.assertIn("ETag", response3)
        self.assertNotEqual(response3["ETag"], etag1, "ETag should change after data modification")

        data3 = response3.json()
        # Verify the modified protein name appears in the response
        spectra = data3["data"]["spectra"]
        # Find any spectrum owned by the modified protein
        names = {s["owner"]["name"] for s in spectra}
        self.assertIn(
            "ModifiedProtein",
            names,
            "Should find spectra with modified protein name",
        )
