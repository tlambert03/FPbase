import pytest
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.urls import reverse
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select

from proteins.models.protein import Protein

SEQ = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
# reverse translation of DGDVNGHKFSVSGEGEGDATYGKLTLKFICT
cDNA = "gatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcacc"


@pytest.mark.usefixtures("uses_frontend", "use_real_webpack_loader")
class TestPagesRender(StaticLiveServerTestCase):
    browser: webdriver.Chrome
    p1: Protein

    @classmethod
    def setUpClass(cls) -> None:
        cls.browser = webdriver.Chrome()
        return super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.browser.quit()
        return super().tearDownClass()

    def setUp(self) -> None:
        self.p1, _ = Protein.objects.get_or_create(name="mySpecialProtein", seq=SEQ)

    def _load_reverse(self, url_name):
        self.browser.get(self.live_server_url + reverse(url_name))

    def test_1spectra(self):
        self._load_reverse("proteins:spectra")
        assert self.browser.get_log("browser") == []

    def test_microscopes(self):
        self._load_reverse("proteins:microscopes")
        assert self.browser.get_log("browser") == []

    def test_blast(self):
        self._load_reverse("proteins:blast")
        assert self.browser.get_log("browser") == []
        text_input = self.browser.find_element(value="queryInput")
        assert text_input.is_displayed()
        text_input.send_keys(SEQ[5:20].replace("LDG", "LG"))
        submit = self.browser.find_element(by="css selector", value='button[type="submit"]')
        submit.click()
        # wait for the results to load
        self.browser.implicitly_wait(1)

        # assert a table is displayed with the first row containing the protein
        r1 = self.browser.find_element(by="xpath", value="//table/tbody/tr[1]/td[1]/a")
        assert r1.text == self.p1.name
        # it should link to the alignment table
        r1.click()

    def test_fret(self):
        self._load_reverse("proteins:fret")
        assert self.browser.get_log("browser") == []

    def test_search(self):
        self._load_reverse("proteins:search")
        assert self.browser.get_log("browser") == []

        elem = self.browser.find_element(value="filter-select-0")
        elem.send_keys("seq")
        assert Select(elem).first_selected_option.text == "Sequence"
        # tab over to action
        elem.send_keys(Keys.TAB)

        elem = self.browser.switch_to.active_element
        elem.send_keys("cdna")
        assert Select(elem).first_selected_option.text == "cDNA could contain"
        elem.send_keys(Keys.TAB)
        elem.send_keys(cDNA)

        # add another filter
        self.browser.find_element(value="add-row-btn").click()

        elem = self.browser.find_element(value="filter-select-1")
        elem.send_keys("name")
        assert Select(elem).first_selected_option.text == "Name or Alias"
        elem.send_keys(Keys.TAB)
        self.browser.switch_to.active_element.send_keys("cont")
        self.browser.switch_to.active_element.send_keys(Keys.TAB)
        self.browser.switch_to.active_element.send_keys(self.p1.name[2:6])

        self.browser.find_element(by="css selector", value='button[type="submit"]').click()
        assert self.browser.current_url == self.live_server_url + self.p1.get_absolute_url()
        assert self.browser.find_element(by="xpath", value="//h1").text == self.p1.name
        shown_seq = self.browser.find_element(by="class name", value="formatted_aminosquence").text
        assert shown_seq.replace(" ", "") == SEQ
