"""Tests for FPbase template tags."""

from __future__ import annotations

import pytest
from django.template import Context, Template


class TestIconTemplateTag:
    """Tests for the icon template tag."""

    def test_icon_basic(self):
        """Test basic icon rendering."""
        t = Template('{% load fpbase_tags %}{% icon "info" %}')
        html = t.render(Context({}))
        assert "fas" in html
        assert "fa-info-circle" in html
        assert "<i " in html
        assert "</i>" in html

    def test_icon_with_class(self):
        """Test icon rendering with additional CSS class."""
        t = Template('{% load fpbase_tags %}{% icon "info" class_="mr-2" %}')
        html = t.render(Context({}))
        assert "mr-2" in html
        assert "fas fa-info-circle" in html

    def test_icon_with_style(self):
        """Test icon rendering with inline style."""
        t = Template('{% load fpbase_tags %}{% icon "info" style="font-size: 0.8rem;" %}')
        html = t.render(Context({}))
        assert 'style="font-size: 0.8rem;"' in html

    def test_icon_with_aria_hidden(self):
        """Test icon rendering with aria-hidden attribute."""
        t = Template('{% load fpbase_tags %}{% icon "info" aria_hidden="true" %}')
        html = t.render(Context({}))
        assert 'aria-hidden="true"' in html

    def test_icon_multiple_attributes(self):
        """Test icon rendering with multiple attributes."""
        t = Template('{% load fpbase_tags %}{% icon "warning" class_="mr-2" style="color: red;" aria_hidden="true" %}')
        html = t.render(Context({}))
        assert "mr-2" in html
        assert 'style="color: red;"' in html
        assert 'aria-hidden="true"' in html
        assert "fas fa-exclamation-circle" in html

    def test_icon_invalid_name(self):
        """Test that invalid icon names raise ValueError."""
        t = Template('{% load fpbase_tags %}{% icon "nonexistent-icon" %}')
        with pytest.raises(ValueError, match="not found in FPbase icon vocabulary"):
            t.render(Context({}))

    def test_icon_xss_protection_class(self):
        """Test that class_ parameter is escaped to prevent XSS."""
        t = Template('{% load fpbase_tags %}{% icon "info" class_=malicious %}')
        html = t.render(Context({"malicious": 'mr-2" onload="alert(1)'}))
        # The malicious code should be escaped
        assert "onload" not in html or "&quot;" in html
        assert "alert" not in html or "&quot;" in html

    def test_icon_xss_protection_style(self):
        """Test that style parameter is escaped to prevent XSS."""
        t = Template('{% load fpbase_tags %}{% icon "info" style=malicious %}')
        html = t.render(Context({"malicious": 'color:red" onload="alert(1)'}))
        # The malicious code should be escaped
        assert "onload" not in html or "&quot;" in html

    def test_icon_xss_protection_custom_attr(self):
        """Test that custom attributes are escaped to prevent XSS."""
        t = Template('{% load fpbase_tags %}{% icon "info" title=malicious %}')
        html = t.render(Context({"malicious": 'Info" onload="alert(1)'}))
        # The malicious code should be escaped
        assert "onload" not in html or "&quot;" in html

    def test_icon_all_semantic_names(self):
        """Test that all defined semantic icon names render without error."""
        icon_names = [
            "info",
            "warning",
            "alert",
            "help",
            "question",
            "close",
            "remove",
            "menu",
            "grid",
            "search",
            "filter",
            "view",
            "settings",
            "edit",
            "delete",
            "trash",
            "undo",
            "check",
            "success",
            "selected",
            "unselected",
            "heart",
            "add",
            "add-item",
            "download",
            "upload",
            "share",
            "share-square",
            "link",
            "external-link",
            "exchange",
            "book",
            "collection",
            "quote",
            "photo",
            "chart",
            "table",
            "flag",
            "clock",
            "spinner",
            "lightbulb",
            "sun",
            "email",
            "wrench",
            "keyboard",
            "google",
            "twitter",
            "orcid",
        ]
        for icon_name in icon_names:
            t = Template(f'{{% load fpbase_tags %}}{{% icon "{icon_name}" %}}')
            html = t.render(Context({}))
            assert "<i " in html, f"Icon '{icon_name}' failed to render"
            assert "</i>" in html, f"Icon '{icon_name}' failed to render"

    def test_icon_fab_style(self):
        """Test that brand icons use fab style."""
        t = Template('{% load fpbase_tags %}{% icon "google" %}')
        html = t.render(Context({}))
        assert "fab" in html
        assert "fa-google" in html

    def test_icon_far_style(self):
        """Test that regular icons use far style."""
        t = Template('{% load fpbase_tags %}{% icon "unselected" %}')
        html = t.render(Context({}))
        assert "far" in html
        assert "fa-square" in html

    def test_icon_spinner_renders_correctly(self):
        """Test that spinner icon uses correct FontAwesome icon name."""
        t = Template('{% load fpbase_tags %}{% icon "spinner" %}')
        html = t.render(Context({}))
        assert "fa-spinner" in html
        # Make sure we're using fa-spinner icon, not fa-spin (which is a CSS class, not an icon)
        assert 'class="fas fa-spinner"' in html

    def test_icon_heart_exists(self):
        """Test that heart icon is available."""
        t = Template('{% load fpbase_tags %}{% icon "heart" %}')
        html = t.render(Context({}))
        assert "fa-heart" in html

    def test_icon_question_exists(self):
        """Test that question icon is available."""
        t = Template('{% load fpbase_tags %}{% icon "question" %}')
        html = t.render(Context({}))
        assert "fa-question-circle" in html

    def test_icon_boolean_attribute(self):
        """Test that boolean attributes are handled correctly."""
        # Test true boolean
        t = Template('{% load fpbase_tags %}{% icon "info" readonly=True %}')
        html = t.render(Context({}))
        assert "readonly" in html
        assert 'readonly="' not in html  # Should be just 'readonly', not 'readonly="True"'

        # Test false boolean
        t = Template('{% load fpbase_tags %}{% icon "info" readonly=False %}')
        html = t.render(Context({}))
        assert "readonly" not in html

    def test_icon_data_attributes(self):
        """Test that data attributes with underscores are converted to hyphens."""
        t = Template('{% load fpbase_tags %}{% icon "info" data_toggle="tooltip" %}')
        html = t.render(Context({}))
        assert 'data-toggle="tooltip"' in html
