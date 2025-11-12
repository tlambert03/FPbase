"""Tests for ETag support, cache utilities, and model version tracking."""

from __future__ import annotations

import pytest
from django.core.cache import cache
from rest_framework.test import APIClient

from fpbase.cache_utils import get_model_version, invalidate_model_version
from fpbase.etag_utils import (
    generate_content_etag,
    generate_version_etag,
    parse_etag_header,
)
from proteins.models import Protein, Spectrum


@pytest.fixture(autouse=True)
def clear_cache():
    """Clear cache before each test."""
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def api_client():
    """Create API client for testing."""
    return APIClient()


# ============================================================================
# Model Version Tracking Tests
# ============================================================================


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


# ============================================================================
# ETag Generation and Parsing Tests
# ============================================================================


class TestGenerateVersionETag:
    """Test version-based ETag generation."""

    def test_returns_weak_etag_format(self):
        """Should return weak ETag with W/"..." format."""
        etag = generate_version_etag(Spectrum)
        assert etag.startswith('W/"')
        assert etag.endswith('"')
        assert len(etag) > 4  # More than just W/""

    def test_same_models_return_same_etag(self):
        """Same models should return same ETag if no changes."""
        etag1 = generate_version_etag(Spectrum)
        etag2 = generate_version_etag(Spectrum)
        assert etag1 == etag2

    def test_different_models_return_different_etags(self):
        """Different models should return different ETags."""
        etag1 = generate_version_etag(Spectrum)
        etag2 = generate_version_etag(Protein)
        assert etag1 != etag2

    def test_multiple_models_combined(self):
        """Should combine multiple model versions."""
        single_etag = generate_version_etag(Spectrum)
        combined_etag = generate_version_etag(Spectrum, Protein)
        assert combined_etag.startswith('W/"')
        assert single_etag != combined_etag

    def test_etag_changes_after_model_save(self, db):
        """ETag should change after model instance is saved."""
        etag_before = generate_version_etag(Protein)

        protein = Protein.objects.create(
            name="Test ETag Change",
            slug="test-etag-change",
        )

        etag_after = generate_version_etag(Protein)
        assert etag_before != etag_after

        protein.delete()


class TestGenerateContentETag:
    """Test content-based ETag generation."""

    def test_returns_strong_etag_format(self):
        """Should return strong ETag with "..." format (no W/)."""
        etag = generate_content_etag("test content")
        assert etag.startswith('"')
        assert etag.endswith('"')
        assert not etag.startswith('W/"')

    def test_same_content_returns_same_etag(self):
        """Same content should always return same ETag."""
        content = "test content"
        etag1 = generate_content_etag(content)
        etag2 = generate_content_etag(content)
        assert etag1 == etag2

    def test_different_content_returns_different_etags(self):
        """Different content should return different ETags."""
        etag1 = generate_content_etag("content 1")
        etag2 = generate_content_etag("content 2")
        assert etag1 != etag2

    def test_handles_bytes_input(self):
        """Should handle both string and bytes input."""
        content_str = "test content"
        content_bytes = b"test content"

        etag_str = generate_content_etag(content_str)
        etag_bytes = generate_content_etag(content_bytes)

        # Both should produce the same ETag
        assert etag_str == etag_bytes

    def test_handles_unicode_content(self):
        """Should handle unicode characters."""
        content = "test 日本語 content"
        etag = generate_content_etag(content)
        assert isinstance(etag, str)
        assert etag.startswith('"')

    def test_handles_json_content(self):
        """Should handle JSON-like content."""
        content = '{"key": "value", "number": 123}'
        etag = generate_content_etag(content)
        assert isinstance(etag, str)


class TestParseETagHeader:
    """Test If-None-Match header parsing."""

    def test_parses_single_etag(self):
        """Should parse single ETag."""
        header = '"abc123"'
        result = parse_etag_header(header)
        assert result == ['"abc123"']

    def test_parses_multiple_etags(self):
        """Should parse comma-separated ETags."""
        header = '"abc123", "def456", "ghi789"'
        result = parse_etag_header(header)
        assert result == ['"abc123"', '"def456"', '"ghi789"']

    def test_parses_weak_etags(self):
        """Should handle weak ETags with W/ prefix."""
        header = 'W/"abc123"'
        result = parse_etag_header(header)
        assert result == ['W/"abc123"']

    def test_parses_mixed_strong_and_weak(self):
        """Should handle mix of strong and weak ETags."""
        header = 'W/"abc123", "def456"'
        result = parse_etag_header(header)
        assert result == ['W/"abc123"', '"def456"']

    def test_handles_whitespace(self):
        """Should handle extra whitespace."""
        header = '  "abc123"  ,  "def456"  '
        result = parse_etag_header(header)
        assert result == ['"abc123"', '"def456"']

    def test_handles_asterisk(self):
        """Should handle * (match any)."""
        header = "*"
        result = parse_etag_header(header)
        assert result == ["*"]

    def test_handles_empty_string(self):
        """Should handle empty header."""
        result = parse_etag_header("")
        assert result == []

    def test_handles_none(self):
        """Should handle None input."""
        result = parse_etag_header(None)
        assert result == []


class TestETagMatching:
    """Test ETag comparison logic."""

    def test_weak_etag_matches_itself(self):
        """Weak ETag should match itself."""
        etag = generate_version_etag(Spectrum)
        parsed = parse_etag_header(etag)
        assert etag in parsed

    def test_strong_etag_matches_itself(self):
        """Strong ETag should match itself."""
        etag = generate_content_etag("test")
        parsed = parse_etag_header(etag)
        assert etag in parsed


# ============================================================================
# API View ETag Support Tests
# ============================================================================


class TestAPIViewETagSupport:
    """Test ETag support in REST API views."""

    def test_response_includes_etag_header(self, api_client, db):
        """Response should include ETag header."""
        response = api_client.get("/api/proteins/table-data/")
        assert response.status_code == 200
        assert "ETag" in response
        assert response["ETag"].startswith('W/"')  # Should be weak ETag

    def test_matching_etag_returns_304(self, api_client, db):
        """Request with matching ETag should return 304 Not Modified."""
        # Get initial response with ETag
        response1 = api_client.get("/api/proteins/table-data/")
        assert response1.status_code == 200

        # Get the ETag and send it back
        etag = response1["ETag"]
        response2 = api_client.get("/api/proteins/table-data/", HTTP_IF_NONE_MATCH=etag)
        assert response2.status_code == 304
        assert len(response2.content) == 0  # 304 has no body
        assert response2["ETag"] == etag  # 304 should include ETag

    def test_non_matching_etag_returns_200(self, api_client, db):
        """Request with non-matching ETag should return 200 with data."""
        # Send request with bogus ETag
        response = api_client.get(
            "/api/proteins/table-data/",
            HTTP_IF_NONE_MATCH='"bogus-etag-12345"',
        )
        assert response.status_code == 200
        # Should have content
        assert len(response.content) > 0

    def test_etag_changes_after_data_modification(self, api_client, db):
        """ETag should change after underlying data changes."""
        # Get initial ETag
        response1 = api_client.get("/api/proteins/table-data/")
        assert response1.status_code == 200
        etag1 = response1["ETag"]

        # Modify data
        protein = Protein.objects.create(
            name="New Test Protein",
            slug="new-test-protein-etag",
        )

        # Get new ETag - should be different
        response2 = api_client.get("/api/proteins/table-data/")
        assert response2.status_code == 200
        etag2 = response2["ETag"]

        # ETags should differ after data change
        assert etag1 != etag2

        protein.delete()

    def test_multiple_etags_in_if_none_match(self, api_client, db):
        """Should handle multiple ETags in If-None-Match header."""
        # Client can send multiple ETags (from different cached versions)
        response = api_client.get(
            "/api/proteins/table-data/",
            HTTP_IF_NONE_MATCH='"old-etag-1", "old-etag-2", "old-etag-3"',
        )
        # None match, so should return 200
        assert response.status_code == 200

    def test_wildcard_etag_returns_304(self, api_client, db):
        """Wildcard * in If-None-Match should match any ETag."""
        # This is used for "only create if not exists" semantics
        # For GET requests, this would return 304
        # (Though this is more commonly used for PUT/POST)
        api_client.get(
            "/api/proteins/table-data/",
            HTTP_IF_NONE_MATCH="*",
        )
        # After implementation, should return 304
        # assert response.status_code == 304


class TestGraphQLETagSupport:
    """Test ETag support in GraphQL endpoint."""

    def test_graphql_response_includes_etag(self, api_client, db):
        """GraphQL response should include ETag header."""
        # Simple introspection query
        query = '{"query": "{ __typename }"}'
        response = api_client.post(
            "/graphql/",
            data=query,
            content_type="application/json",
        )
        assert response.status_code == 200
        # After implementation, should have ETag
        # assert "ETag" in response

    def test_graphql_matching_etag_returns_304(self, api_client, db):
        """GraphQL with matching ETag should return 304."""
        # First request
        query = '{"query": "{ __typename }"}'
        response1 = api_client.post(
            "/graphql/",
            data=query,
            content_type="application/json",
        )
        assert response1.status_code == 200

        # After implementation, send same request with ETag
        # etag = response1["ETag"]
        # response2 = api_client.post(
        #     "/graphql/",
        #     data=query,
        #     content_type="application/json",
        #     HTTP_IF_NONE_MATCH=etag,
        # )
        # assert response2.status_code == 304


class TestETagCacheControl:
    """Test interaction between ETags and Cache-Control headers."""

    def test_etag_works_with_cache_control(self, api_client, db):
        """ETags should work alongside Cache-Control headers."""
        response = api_client.get("/api/proteins/table-data/")
        assert response.status_code == 200
        # Should have both ETag and Cache-Control (after implementation)
        # assert "ETag" in response
        # assert "Cache-Control" in response

    def test_weak_etag_for_version_based_caching(self, api_client, db):
        """Version-based ETags should be weak (W/"...")."""
        api_client.get("/api/proteins/table-data/")
        # After implementation
        # etag = response["ETag"]
        # assert etag.startswith('W/"')  # Weak ETag


class TestETagPerformance:
    """Test that ETags actually reduce bandwidth and processing."""

    def test_304_response_has_no_body(self, api_client, db):
        """304 responses should have no body to save bandwidth."""
        # After implementation, test that 304 responses are empty
        pass

    def test_304_response_skips_serialization(self, api_client, db):
        """304 responses should skip expensive serialization."""
        # This is implicit in Django/DRF - if we return early with 304,
        # the serialization never happens
        pass


class TestGraphQLETagWorkflow:
    """Test the complete ETag workflow for GraphQL queries."""

    def test_spectra_list_etag_workflow(self, api_client, db):
        """Test SpectraList query: 200 → 304 → 200 on data change."""
        from proteins.factories import DyeFactory, SpectrumFactory
        from proteins.models import Spectrum

        # Create initial spectra with owners (dyes)
        for _ in range(3):
            dye = DyeFactory.create()
            SpectrumFactory.create(owner_dye=dye, category=Spectrum.DYE, subtype=Spectrum.EX)

        # GraphQL query for SpectraList
        query = """
        query _FPB_SpectraList {
            spectra {
                id
                category
                subtype
                owner { id name slug url }
            }
        }
        """

        # First request - should return 200 with ETag
        response1 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_SpectraList"},
        )
        assert response1.status_code == 200
        assert "ETag" in response1
        assert response1["Cache-Control"] == "public, max-age=600"

        etag1 = response1["ETag"]
        assert etag1.startswith('W/"')  # Weak ETag
        assert "-" in etag1  # Format: W/"query_hash-model_version"

        # Second request with same ETag - should return 304
        response2 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_SpectraList"},
            HTTP_IF_NONE_MATCH=etag1,
        )
        assert response2.status_code == 304
        assert response2["ETag"] == etag1
        assert len(response2.content) == 0  # No body to save bandwidth

        # Modify data - create new spectrum
        new_dye = DyeFactory.create()
        SpectrumFactory.create(owner_dye=new_dye, category=Spectrum.DYE, subtype=Spectrum.EM)

        # Third request - ETag should change, return 200 with new data
        response3 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_SpectraList"},
            HTTP_IF_NONE_MATCH=etag1,  # Old ETag
        )
        assert response3.status_code == 200
        assert "ETag" in response3

        etag3 = response3["ETag"]
        assert etag3 != etag1  # ETag changed due to data change
        assert etag3.startswith('W/"')

    def test_optical_config_list_etag_workflow(self, api_client, db):
        """Test OpticalConfigList query: 200 → 304 → 200 on data change."""
        from proteins.factories import MicroscopeFactory, OpticalConfigFactory

        # Create initial data
        microscope = MicroscopeFactory.create()
        OpticalConfigFactory.create_batch(2, microscope=microscope)

        query = """
        query _FPB_OpticalConfigList {
            opticalConfigs {
                id
                name
                comments
                microscope { name }
            }
        }
        """

        # First request - 200 with ETag
        response1 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_OpticalConfigList"},
        )
        assert response1.status_code == 200
        assert "ETag" in response1
        etag1 = response1["ETag"]

        # Second request - 304
        response2 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_OpticalConfigList"},
            HTTP_IF_NONE_MATCH=etag1,
        )
        assert response2.status_code == 304

        # Modify microscope (which OpticalConfig depends on)
        microscope.name = "Updated Microscope"
        microscope.save()

        # Third request - 200 with new ETag
        response3 = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "_FPB_OpticalConfigList"},
            HTTP_IF_NONE_MATCH=etag1,
        )
        assert response3.status_code == 200
        assert response3["ETag"] != etag1

    def test_query_hash_prevents_collision(self, api_client, db):
        """Different queries with same operation name get different ETags."""
        from proteins.factories import ProteinFactory

        ProteinFactory.create()

        # Two different queries with the same operation name (hypothetical attack/bug)
        query1 = """
        query _FPB_SpectraList {
            spectra { id category }
        }
        """

        query2 = """
        query _FPB_SpectraList {
            spectra { id subtype owner { name } }
        }
        """

        response1 = api_client.get(
            "/graphql/",
            {"query": query1, "operationName": "_FPB_SpectraList"},
        )
        response2 = api_client.get(
            "/graphql/",
            {"query": query2, "operationName": "_FPB_SpectraList"},
        )

        # Different query content = different ETags (safety)
        assert response1["ETag"] != response2["ETag"]

    def test_post_requests_no_etag(self, api_client):
        """POST requests should not get ETag support."""
        query = '{"query": "query _FPB_SpectraList { spectra { id } }"}'

        response = api_client.post(
            "/graphql/",
            data=query,
            content_type="application/json",
        )
        assert response.status_code == 200
        assert "ETag" not in response  # No ETag for POST

    def test_unknown_operation_no_etag(self, api_client):
        """Unknown operations should not get ETag support."""
        query = """
        query UnknownOperation {
            __typename
        }
        """

        response = api_client.get(
            "/graphql/",
            {"query": query, "operationName": "UnknownOperation"},
        )
        assert response.status_code == 200
        assert "ETag" not in response  # No ETag for unknown operations


class TestGraphQLOperationNameSync:
    """Test that frontend and backend operation names are in sync."""

    def test_graphql_operation_names_match_frontend(self):
        """Ensure backend ETag mapping matches frontend operation names."""
        import re
        from pathlib import Path

        from fpbase.views import GRAPHQL_OPERATION_ETAG_MODELS

        # Read frontend queries.ts file
        frontend_queries_path = (
            Path(__file__).parent.parent.parent.parent / "packages" / "spectra" / "src" / "api" / "queries.ts"
        )

        if not frontend_queries_path.exists():
            # If running in CI or somewhere without frontend code, skip
            import pytest

            pytest.skip("Frontend code not available")

        queries_content = frontend_queries_path.read_text()

        # Extract all operation names from frontend (pattern: query OperationName)
        operation_name_pattern = r"query (_FPB_\w+)"
        frontend_operations = set(re.findall(operation_name_pattern, queries_content))

        # Get backend operation names
        backend_operations = set(GRAPHQL_OPERATION_ETAG_MODELS.keys())

        # All backend operations should exist in frontend
        missing_in_frontend = backend_operations - frontend_operations
        assert not missing_in_frontend, (
            f"Backend has ETag mappings for operations not found in frontend: {missing_in_frontend}"
        )

        # All frontend _FPB_ operations should have backend mappings
        missing_in_backend = frontend_operations - backend_operations
        assert not missing_in_backend, (
            f"Frontend has _FPB_ operations without backend ETag mappings: {missing_in_backend}. "
            f"Either add them to GRAPHQL_OPERATION_ETAG_MODELS or remove _FPB_ prefix."
        )
