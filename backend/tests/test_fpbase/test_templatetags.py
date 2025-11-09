"""Tests for FPbase template tags."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pytest
from django.template import Context, Template

from fpbase.templatetags.fpbase_tags import AVAILABLE_ICONS


class TestIconTemplateTag:
    """Tests for the icon template tag."""

    def test_icon_basic(self):
        """Test basic icon rendering."""
        t = Template('{% load fpbase_tags %}{% icon "info" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "</svg>" in html
        assert "viewBox=" in html
        assert "<path" in html

    def test_icon_with_class(self):
        """Test icon rendering with additional CSS class."""
        t = Template('{% load fpbase_tags %}{% icon "info" class_="mr-2" %}')
        html = t.render(Context({}))
        assert 'class="mr-2"' in html
        assert "<svg" in html

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
        assert 'class="mr-2"' in html
        assert 'style="color: red;"' in html
        assert 'aria-hidden="true"' in html
        assert "<svg" in html

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
        for icon_name in AVAILABLE_ICONS:
            t = Template(f'{{% load fpbase_tags %}}{{% icon "{icon_name}" %}}')
            html = t.render(Context({}))
            assert "<svg" in html, f"Icon '{icon_name}' failed to render"
            assert "</svg>" in html, f"Icon '{icon_name}' failed to render"
            assert "viewBox=" in html, f"Icon '{icon_name}' missing viewBox"
            assert "<path" in html, f"Icon '{icon_name}' missing path"

    def test_icon_brand_renders(self):
        """Test that brand icons render as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "google" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_unselected_renders(self):
        """Test that unselected icon renders as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "unselected" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_spinner_renders_correctly(self):
        """Test that spinner icon renders as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "spinner" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_heart_exists(self):
        """Test that heart icon is available."""
        t = Template('{% load fpbase_tags %}{% icon "heart" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_question_exists(self):
        """Test that question icon is available."""
        t = Template('{% load fpbase_tags %}{% icon "question" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

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

    def test_icon_keyboard_renders(self):
        """Test that keyboard icon renders as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "keyboard" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_flag_outline_renders(self):
        """Test that flag-outline icon renders as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "flag-outline" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_icon_flag_solid_renders(self):
        """Test that flag icon renders as SVG."""
        t = Template('{% load fpbase_tags %}{% icon "flag" %}')
        html = t.render(Context({}))
        assert "<svg" in html
        assert "<path" in html

    def test_all_template_icons_are_valid(self):
        """Test that all icon names used in templates exist in AVAILABLE_ICONS.

        This test scans all .html templates in the backend directory and verifies
        that every {% icon "..." %} tag uses a valid icon name from AVAILABLE_ICONS.
        """
        # Get the backend directory (parent of tests directory)
        backend_dir = Path(__file__).parent.parent.parent

        # Pattern to match {% icon "name" %} or {% icon 'name' %}
        # This matches both double and single quotes
        icon_pattern = re.compile(r'{%\s*icon\s+["\']([^"\']+)["\']')

        invalid_icons = []
        icon_usages = defaultdict(list)

        # Find all .html files
        for template_file in backend_dir.rglob("*.html"):
            content = template_file.read_text(encoding="utf-8")

            # Find all icon tags in this file
            for match in icon_pattern.finditer(content):
                icon_name = match.group(1)
                icon_usages[icon_name].append(str(template_file.relative_to(backend_dir)))

                # Check if icon exists in AVAILABLE_ICONS
                if icon_name not in AVAILABLE_ICONS:
                    # Get line number for better error reporting
                    line_num = content[: match.start()].count("\n") + 1
                    invalid_icons.append(
                        {
                            "icon": icon_name,
                            "file": str(template_file.relative_to(backend_dir)),
                            "line": line_num,
                        }
                    )

        # If there are invalid icons, create a helpful error message
        if invalid_icons:
            error_lines = ["Found icon tags with invalid icon names:"]
            for item in invalid_icons:
                error_lines.append(f"  - '{item['icon']}' in {item['file']}:{item['line']}")
            error_lines.append("")
            error_lines.append("Available icons in AVAILABLE_ICONS:")
            error_lines.append(f"  {', '.join(sorted(AVAILABLE_ICONS))}")

            pytest.fail("\n".join(error_lines))

        # Also verify we found at least some icon usages (sanity check)
        assert len(icon_usages) > 0, (
            "No {% icon %} tags found in templates. This test may not be searching the correct directory."
        )
