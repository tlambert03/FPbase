"""Tests for ETag generation and parsing utilities."""

from __future__ import annotations

import pytest
from django.core.cache import cache

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
