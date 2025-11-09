"""Tests for model version tracking and cache utilities."""

from __future__ import annotations

import pytest
from django.core.cache import cache

from fpbase.cache_utils import get_model_version, invalidate_model_version
from proteins.models import Protein, Spectrum


@pytest.fixture(autouse=True)
def clear_cache():
    """Clear cache before each test."""
    cache.clear()
    yield
    cache.clear()


class TestModelVersionTracking:
    """Test model version tracking system."""

    def test_get_model_version_returns_string(self):
        """get_model_version should return a string hash."""
        version = get_model_version(Spectrum)
        assert isinstance(version, str)
        assert len(version) > 0

    def test_get_model_version_is_consistent(self):
        """Same models should return same version if no changes."""
        version1 = get_model_version(Spectrum)
        version2 = get_model_version(Spectrum)
        assert version1 == version2

    def test_get_model_version_combines_multiple_models(self):
        """Should combine versions from multiple models."""
        version = get_model_version(Spectrum, Protein)
        assert isinstance(version, str)
        # Combined version should be different from single model
        assert version != get_model_version(Spectrum)

    def test_invalidate_model_version_changes_version(self):
        """Invalidating should change the version."""
        version_before = get_model_version(Spectrum)
        invalidate_model_version(Spectrum)
        version_after = get_model_version(Spectrum)
        assert version_before != version_after

    def test_invalidate_only_affects_specified_model(self):
        """Invalidating one model shouldn't affect others."""
        spectrum_v1 = get_model_version(Spectrum)
        protein_v1 = get_model_version(Protein)

        invalidate_model_version(Spectrum)

        spectrum_v2 = get_model_version(Spectrum)
        protein_v2 = get_model_version(Protein)

        assert spectrum_v1 != spectrum_v2  # Changed
        assert protein_v1 == protein_v2  # Unchanged

    def test_model_save_invalidates_version(self, db):
        """Saving a model instance should automatically invalidate version."""
        version_before = get_model_version(Protein)

        # Create a new protein (simpler model with fewer validation requirements)
        protein = Protein.objects.create(
            name="Test Protein",
            slug="test-protein-etag",
        )

        version_after = get_model_version(Protein)
        assert version_before != version_after

        # Cleanup
        protein.delete()

    def test_model_delete_invalidates_version(self, db):
        """Deleting a model instance should automatically invalidate version."""
        # Create a protein
        protein = Protein.objects.create(
            name="Test Protein Delete",
            slug="test-protein-delete-etag",
        )

        version_before = get_model_version(Protein)
        protein.delete()
        version_after = get_model_version(Protein)

        assert version_before != version_after

    def test_cache_key_format(self):
        """Cache keys should use app_label.ModelName format."""
        # This is more of an implementation detail test
        # We'll access the internal key format
        get_model_version(Spectrum)
        # The key should be stored in cache
        cached = cache.get("model_version:proteins.Spectrum")
        assert cached is not None


class TestVersionCaching:
    """Test that versions are properly cached."""

    def test_version_persists_in_cache(self):
        """Version should persist in cache between calls."""
        get_model_version(Spectrum)  # First call, sets cache
        cache_key = "model_version:proteins.Spectrum"
        cached_value = cache.get(cache_key)
        assert cached_value is not None

    def test_multiple_model_version_uses_all_caches(self):
        """get_model_version with multiple models should read all caches."""
        # Set up known versions
        get_model_version(Spectrum)
        get_model_version(Protein)

        # Both should be cached
        assert cache.get("model_version:proteins.Spectrum") is not None
        assert cache.get("model_version:proteins.Protein") is not None

        # Combined version should use both
        combined = get_model_version(Spectrum, Protein)
        assert isinstance(combined, str)
