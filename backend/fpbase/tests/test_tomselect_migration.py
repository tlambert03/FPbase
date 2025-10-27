"""Test stubs for Tom-Select migration following TDD approach."""

from __future__ import annotations

import pytest
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium import webdriver
from selenium.webdriver.common.by import By


@pytest.mark.usefixtures("uses_frontend", "use_real_webpack_loader")
class TestTomSelectMigration(StaticLiveServerTestCase):
    """Tests for Tom-Select migration progress."""

    browser: webdriver.Chrome

    @classmethod
    def setUpClass(cls) -> None:
        options = webdriver.ChromeOptions()
        cls.browser = webdriver.Chrome(options=options)
        return super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.browser.quit()
        return super().tearDownClass()

    @pytest.mark.xfail(reason="django-tomselect URL resolution causes circular import during test setup")
    def test_microscope_form_uses_tomselect(self):
        """Test that microscope form uses Tom-Select instead of Select2."""
        self.browser.get(f"{self.live_server_url}/microscope/create/")
        # Check that Tom-Select wrapper exists
        ts_wrappers = self.browser.find_elements(By.CLASS_NAME, "ts-wrapper")
        assert len(ts_wrappers) > 0, "Tom-Select wrapper not found"
        # Check that Select2 container does not exist
        select2_containers = self.browser.find_elements(By.CLASS_NAME, "select2-container")
        assert len(select2_containers) == 0, "Select2 container still exists"

    def test_lineage_form_uses_tomselect_field(self):
        """Test that LineageForm uses TomSelectModelChoiceField."""
        from django_tomselect.forms import TomSelectModelChoiceField

        from proteins.forms.forms import LineageForm

        form = LineageForm()
        assert isinstance(form.fields["parent"], TomSelectModelChoiceField)
