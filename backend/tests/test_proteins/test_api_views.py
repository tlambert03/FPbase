"""
Tests for protein API views.

These tests ensure the API endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import MicroscopeFactory, OpticalConfigFactory, ProteinFactory, StateFactory


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


class OpticalConfigListAPIViewTests(TestCase):
    """Test the optical_configs_list endpoint with ETag support."""

    @classmethod
    def setUpTestData(cls):
        """Create test optical configs."""
        microscope = MicroscopeFactory(name="Test Microscope")
        for i in range(3):
            OpticalConfigFactory(microscope=microscope, name=f"Config {i}")

    def test_optical_configs_list_returns_data(self):
        """Test that the endpoint returns optical configs data."""
        response = self.client.get("/api/optical-configs-list/")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "application/json")
        data = response.json()

        self.assertIn("data", data)
        self.assertIn("opticalConfigs", data["data"])
        self.assertEqual(len(data["data"]["opticalConfigs"]), 3)

        # Check structure of returned data
        config = data["data"]["opticalConfigs"][0]
        self.assertIn("id", config)
        self.assertIn("name", config)
        self.assertIn("comments", config)
        self.assertIn("microscope", config)
        self.assertIn("id", config["microscope"])
        self.assertIn("name", config["microscope"])

    def test_optical_configs_list_etag_header(self):
        """Test that the endpoint returns an ETag header."""
        response = self.client.get("/api/optical-configs-list/")

        self.assertEqual(response.status_code, 200)
        self.assertIn("ETag", response)
        self.assertTrue(response["ETag"].startswith('W/"'))
        self.assertIn("Cache-Control", response)
        self.assertIn("max-age=0", response["Cache-Control"])
        self.assertIn("must-revalidate", response["Cache-Control"])

    def test_optical_configs_list_etag_304(self):
        """Test that the endpoint returns 304 when ETag matches."""
        # First request to get the ETag
        response1 = self.client.get("/api/optical-configs-list/")
        etag = response1["ETag"]

        # Second request with If-None-Match header
        response2 = self.client.get("/api/optical-configs-list/", headers={"If-None-Match": etag})

        self.assertEqual(response2.status_code, 304)
        self.assertEqual(response2.content, b"")
        self.assertEqual(response2["ETag"], etag)
