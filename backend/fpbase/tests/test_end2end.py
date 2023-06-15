import pytest
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.urls import reverse
from selenium import webdriver


@pytest.mark.usefixtures("uses_frontend", "use_real_webpack_loader")
class TestPagesRender(StaticLiveServerTestCase):
    def setUp(self):
        self.browser = webdriver.Chrome()

    def _load_reverse(self, url_name):
        self.browser.get(self.live_server_url + reverse(url_name))

    def test_spectra(self):
        self._load_reverse("proteins:spectra")
        assert self.browser.get_log("browser") == []

    def test_microscopes(self):
        self._load_reverse("proteins:microscopes")
        assert self.browser.get_log("browser") == []

    def test_blast(self):
        self._load_reverse("proteins:blast")
        assert self.browser.get_log("browser") == []

    def test_fret(self):
        self._load_reverse("proteins:fret")
        assert self.browser.get_log("browser") == []

    def test_search(self):
        self._load_reverse("proteins:search")
        assert self.browser.get_log("browser") == []
