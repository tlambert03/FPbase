"""
Tests for ajax views in proteins app.

These tests ensure ajax endpoints perform efficiently and avoid N+1 query issues.
"""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from proteins.factories import (
    CameraFactory,
    DyeFactory,
    FilterFactory,
    LightFactory,
    ProteinFactory,
    SpectrumFactory,
    StateFactory,
)
from proteins.models import Spectrum
from proteins.models.dye import DyeState

User = get_user_model()


class SimilarSpectrumOwnersViewTests(TestCase):
    """Test the similar_spectrum_owners ajax endpoint."""

    @classmethod
    def setUpTestData(cls):
        """Create test data once for all tests in this class."""
        # Create proteins with states - all with similar names so find_similar_owners finds them
        cls.proteins = []
        cls.states = []
        for i in range(5):
            protein = ProteinFactory(name=f"SimilarProtein{i}")
            state = StateFactory(protein=protein, name=f"state{i}")
            cls.proteins.append(protein)
            cls.states.append(state)
            # Add spectra to states
            SpectrumFactory(owner_fluor=state, category=Spectrum.PROTEIN, subtype=Spectrum.EX)
            SpectrumFactory(owner_fluor=state, category=Spectrum.PROTEIN, subtype=Spectrum.EM)

        # Create dyes - note: DyeState is the fluorophore that owns spectra, not Dye itself
        cls.dyes = []
        cls.dye_states = []
        for i in range(3):
            dye = DyeFactory(name=f"SimilarDye{i}")
            cls.dyes.append(dye)
            # Create DyeState as the actual fluorophore owner of spectra
            dye_state = DyeState.objects.create(
                dye=dye,
                name=f"SimilarDye{i} state",
                slug=f"similardye{i}-state",
                ex_max=488,
                em_max=520,
            )
            cls.dye_states.append(dye_state)
            # Add spectra to dye states
            SpectrumFactory(owner_fluor=dye_state, category=Spectrum.DYE, subtype=Spectrum.EX)
            SpectrumFactory(owner_fluor=dye_state, category=Spectrum.DYE, subtype=Spectrum.EM)

        # Create filters (factories automatically create spectrum)
        cls.filters = []
        for i in range(2):
            filt = FilterFactory(name=f"SimilarFilter{i}")
            cls.filters.append(filt)

        # Create lights (factories automatically create spectrum)
        cls.lights = []
        for i in range(2):
            light = LightFactory(name=f"SimilarLight{i}")
            cls.lights.append(light)

        # Create cameras (factories automatically create spectrum)
        cls.cameras = []
        for i in range(2):
            camera = CameraFactory(name=f"SimilarCamera{i}")
            cls.cameras.append(camera)

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_states(self):
        """
        Test that similar_spectrum_owners doesn't generate N+1 queries when returning States.

        With proper prefetching, the query count should be constant regardless of the
        number of similar owners found.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "SimilarProtein"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find proteins (limited to 4)
        self.assertGreater(len(data["similars"]), 0)
        self.assertLessEqual(len(data["similars"]), 4)

        query_count = len(context.captured_queries)
        # Expected queries:
        # find_similar_owners makes ~5 queries (one per owner type)
        # Then we fetch each type with proper prefetching (~3-5 queries)
        # Total should be around 20-25 queries
        self.assertLessEqual(
            query_count,
            30,
            f"similar_spectrum_owners generated {query_count} queries for States. "
            f"Expected ≤30 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_dyes(self):
        """
        Test that similar_spectrum_owners doesn't generate N+1 queries when returning Dyes.

        With proper prefetching, the query count should be constant.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "SimilarDye"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find dyes
        self.assertGreater(len(data["similars"]), 0)

        query_count = len(context.captured_queries)
        # Expected queries similar to states test
        self.assertLessEqual(
            query_count,
            30,
            f"similar_spectrum_owners generated {query_count} queries for Dyes. "
            f"Expected ≤30 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_handles_filters(self):
        """Test that Filters (which have 'spectrum' not 'spectra') work correctly."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "SimilarFilter"},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find filters
        self.assertGreater(len(data["similars"]), 0)

        # Check data structure
        for similar in data["similars"]:
            self.assertIn("slug", similar)
            self.assertIn("name", similar)
            self.assertIn("url", similar)
            self.assertIn("spectra", similar)
            self.assertIsInstance(similar["spectra"], list)
            # Filter should have at least one spectrum type
            self.assertGreater(len(similar["spectra"]), 0)

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_handles_lights(self):
        """Test that Lights (which have 'spectrum' not 'spectra') work correctly."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "SimilarLight"},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find lights
        self.assertGreater(len(data["similars"]), 0)

        # Check data structure
        for similar in data["similars"]:
            self.assertIn("spectra", similar)
            self.assertIsInstance(similar["spectra"], list)

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_handles_cameras(self):
        """Test that Cameras (which have 'spectrum' not 'spectra') work correctly."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "SimilarCamera"},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find cameras
        self.assertGreater(len(data["similars"]), 0)

        # Check data structure
        for similar in data["similars"]:
            self.assertIn("spectra", similar)
            self.assertIsInstance(similar["spectra"], list)

    @override_settings(DEBUG=True)
    def test_similar_spectrum_owners_query_count_with_mixed_types(self):
        """
        Test query efficiency with mixed types including Filters/Lights/Cameras.

        This tests the most complex scenario where different types need different prefetching.
        """
        with CaptureQueriesContext(connection) as context:
            response = self.client.post(
                "/ajax/validate_spectrumownername/",
                {"owner": "Similar"},
                headers={"X-Requested-With": "XMLHttpRequest"},
            )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        # Should find multiple items (limited to 4)
        self.assertGreater(len(data["similars"]), 0)
        self.assertLessEqual(len(data["similars"]), 4)

        query_count = len(context.captured_queries)
        # Expected queries:
        # find_similar_owners queries for each owner type (5 types)
        # Then fetch queries for each type found with proper prefetching
        self.assertLessEqual(
            query_count,
            30,
            f"similar_spectrum_owners generated {query_count} queries for mixed types. "
            f"Expected ≤30 queries with proper prefetch_related. "
            f"This may indicate an N+1 query regression.",
        )

    def test_similar_spectrum_owners_returns_correct_data_structure(self):
        """Test that the endpoint returns the expected data structure."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "SimilarProtein"},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        self.assertEqual(response.status_code, 200)
        data = response.json()

        # Check structure
        self.assertIn("similars", data)
        self.assertGreater(len(data["similars"]), 0)

        # Check each similar has expected fields
        for similar in data["similars"]:
            self.assertIn("slug", similar)
            self.assertIn("name", similar)
            self.assertIn("url", similar)
            self.assertIn("spectra", similar)
            self.assertIsInstance(similar["spectra"], list)

    def test_similar_spectrum_owners_state_includes_protein_name(self):
        """Test that States return their protein's name, not the state name."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": self.proteins[0].name},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        data = response.json()
        self.assertGreater(len(data["similars"]), 0)
        # At least one should be the protein we searched for
        names = [s["name"] for s in data["similars"]]
        self.assertIn(self.proteins[0].name, names)

    def test_similar_spectrum_owners_dye_uses_own_name(self):
        """Test that DyeStates (Fluorophores without protein) use their label."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": self.dye_states[0].label},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        data = response.json()
        self.assertGreater(len(data["similars"]), 0)
        # At least one should be the dye state we searched for
        names = [s["name"] for s in data["similars"]]
        self.assertIn(self.dye_states[0].label, names)

    def test_similar_spectrum_owners_works_without_ajax_header(self):
        """Test that the endpoint works without the X-Requested-With header."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "Test"},
        )

        # Should work fine without AJAX header (JSON-only endpoint)
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("similars", data)

    def test_similar_spectrum_owners_limits_to_four_results(self):
        """Test that the endpoint limits results to 4 items."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "Similar"},  # Should match many items
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        data = response.json()
        # Should only return 4 even though more are available
        self.assertLessEqual(len(data["similars"]), 4)

    def test_similar_spectrum_owners_empty_results(self):
        """Test handling of empty results from find_similar_owners."""
        response = self.client.post(
            "/ajax/validate_spectrumownername/",
            {"owner": "NonExistentOwnerNameThatWillNeverMatch12345"},
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(len(data["similars"]), 0)
