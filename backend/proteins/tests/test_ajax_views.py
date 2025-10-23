"""
Tests for ajax views in proteins app.

These tests ensure ajax endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import DyeFactory, ProteinFactory, SpectrumFactory, StateFactory
from proteins.models import Spectrum

User = get_user_model()


class SimilarSpectrumOwnersViewTests(TestCase):
    """Test the similar_spectrum_owners ajax endpoint."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests in this class."""
        # Create proteins with states
        cls.proteins = []
        cls.states = []
        for i in range(5):
            protein = ProteinFactory(name=f"TestProtein{i}")
            state = StateFactory(protein=protein, name=f"state{i}")
            cls.proteins.append(protein)
            cls.states.append(state)
            # Add spectra to states
            for _ in range(2):
                SpectrumFactory(owner_state=state, category=Spectrum.PROTEIN, subtype=Spectrum.EX)

        # Create dyes (Fluorophores without protein)
        cls.dyes = []
        for i in range(3):
            dye = DyeFactory(name=f"TestDye{i}")
            cls.dyes.append(dye)
            # Add spectra to dyes
            for _ in range(2):
                SpectrumFactory(owner_dye=dye, category=Spectrum.DYE, subtype=Spectrum.EX)

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_states(self):
        """
        Test that similar_spectrum_owners doesn't generate N+1 queries when returning States.

        With proper prefetching, the query count should be constant regardless of the
        number of similar owners found.
        """
        # Mock find_similar_owners to return States
        mock_similars = self.states[:4]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            with CaptureQueriesContext(connection) as context:
                response = self.client.post(
                    "/ajax/validate_spectrumownername/",
                    {"owner": "TestProtein"},
                    headers={"X-Requested-With": "XMLHttpRequest"},
                )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data["similars"]), 4)

        query_count = len(context.captured_queries)
        # Expected queries:
        # 1. find_similar_owners (mocked but still executes)
        # 2. Fetch States with protein
        # 3. Prefetch spectra for States
        self.assertLessEqual(
            query_count,
            5,
            f"similar_spectrum_owners generated {query_count} queries for States. "
            f"Expected ≤5 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_dyes(self):
        """
        Test that similar_spectrum_owners doesn't generate N+1 queries when returning Dyes.

        With proper prefetching, the query count should be constant.
        """
        # Mock find_similar_owners to return Dyes (Fluorophores)
        mock_similars = self.dyes[:3]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            with CaptureQueriesContext(connection) as context:
                response = self.client.post(
                    "/ajax/validate_spectrumownername/",
                    {"owner": "TestDye"},
                    headers={"X-Requested-With": "XMLHttpRequest"},
                )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data["similars"]), 3)

        query_count = len(context.captured_queries)
        # Expected queries:
        # 1. find_similar_owners (mocked)
        # 2. Fetch Fluorophores with spectra
        self.assertLessEqual(
            query_count,
            4,
            f"similar_spectrum_owners generated {query_count} queries for Dyes. "
            f"Expected ≤4 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_mixed_types(self):
        """
        Test query efficiency with mixed States and Fluorophores.

        This tests the most complex scenario where different types need different prefetching.
        """
        # Mock find_similar_owners to return a mix
        mock_similars = [self.states[0], self.dyes[0], self.states[1], self.dyes[1]]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            with CaptureQueriesContext(connection) as context:
                response = self.client.post(
                    "/ajax/validate_spectrumownername/",
                    {"owner": "Test"},
                    headers={"X-Requested-With": "XMLHttpRequest"},
                )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data["similars"]), 4)

        query_count = len(context.captured_queries)
        # Expected queries:
        # 1. find_similar_owners (mocked)
        # 2. Fetch States with protein
        # 3. Prefetch spectra for States
        # 4. Fetch Fluorophores (Dyes) with spectra
        self.assertLessEqual(
            query_count,
            6,
            f"similar_spectrum_owners generated {query_count} queries for mixed types. "
            f"Expected ≤6 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    def test_similar_spectrum_owners_returns_correct_data_structure(self):
        """Test that the endpoint returns the expected data structure."""
        # Mock find_similar_owners to return States
        mock_similars = self.states[:2]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "TestProtein"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        self.assertEqual(response.status_code, 200)
        data = response.json()

        # Check structure
        self.assertIn("similars", data)
        self.assertEqual(len(data["similars"]), 2)

        # Check each similar has expected fields
        for similar in data["similars"]:
            self.assertIn("slug", similar)
            self.assertIn("name", similar)
            self.assertIn("url", similar)
            self.assertIn("spectra", similar)
            self.assertIsInstance(similar["spectra"], list)

    def test_similar_spectrum_owners_state_includes_protein_name(self):
        """Test that States return their protein's name, not the state name."""
        state = self.states[0]
        mock_similars = [state]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "TestProtein"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        data = response.json()
        self.assertEqual(data["similars"][0]["name"], state.protein.name)

    def test_similar_spectrum_owners_dye_uses_own_name(self):
        """Test that Dyes (Fluorophores without protein) use their own name."""
        dye = self.dyes[0]
        mock_similars = [dye]

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "TestDye"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        data = response.json()
        self.assertEqual(data["similars"][0]["name"], dye.name)

    def test_similar_spectrum_owners_requires_ajax(self):
        """Test that the endpoint requires an AJAX request."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "Test"},
        )

        # Should return 405 Method Not Allowed for non-AJAX requests
        self.assertEqual(response.status_code, 405)

    def test_similar_spectrum_owners_limits_to_four_results(self):
        """Test that the endpoint limits results to 4 items."""
        # Create more than 4 similar items
        mock_similars = self.states[:5]  # 5 items

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "TestProtein"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        data = response.json()
        # Should only return 4 even though 5 were available
        self.assertEqual(len(data["similars"]), 4)

    def test_similar_spectrum_owners_empty_results(self):
        """Test handling of empty results from find_similar_owners."""
        mock_similars = []

        with patch.object(Spectrum.objects, "find_similar_owners", return_value=mock_similars):
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "NonExistent"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data["similars"]), 0)
