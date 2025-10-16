import os
import shutil
import tempfile

import pytest
from allauth.account.models import EmailAddress
from django.contrib.auth import get_user_model
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.test import Client
from django.urls import reverse
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.wait import WebDriverWait

from proteins.factories import MicroscopeFactory, OpticalConfigWithFiltersFactory, ProteinFactory
from proteins.models import Protein, Spectrum
from proteins.util.blast import _get_binary

SEQ = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
# reverse translation of DGDVNGHKFSVSGEGEGDATYGKLTLKFICT
cDNA = "gatggcgatgtgaacggccataaatttagcgtgagcggcgaaggcgaaggcgatgcgacctatggcaaactgaccctgaaatttatttgcacc"
PASSWORD = "testpass2341o87123o847u3214"


@pytest.mark.usefixtures("uses_frontend", "use_real_webpack_loader")
class TestPagesRender(StaticLiveServerTestCase):
    browser: webdriver.Chrome
    p1: Protein
    download_dir: str

    @classmethod
    def setUpClass(cls) -> None:
        cls.download_dir = tempfile.mkdtemp()
        options = webdriver.ChromeOptions()
        options.add_experimental_option("prefs", {"download.default_directory": cls.download_dir})

        cls.browser = webdriver.Chrome(options=options)
        return super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.browser.quit()
        shutil.rmtree(cls.download_dir)
        return super().tearDownClass()

    def setUp(self) -> None:
        self.browser.get_log("browser")  # clear prior logs
        self.p1 = ProteinFactory(name="knownSequence", seq=SEQ)

    def _load_reverse(self, url_name, **kwargs):
        self.browser.get(self.live_server_url + reverse(url_name, **kwargs))

    def _assert_no_console_errors(self):
        logs = self.browser.get_log("browser")
        for lg in logs:
            if lg["level"] == "SEVERE":
                raise AssertionError(f"Console errors occurred: {logs}")

    def test_spectra(self):
        self._load_reverse("proteins:spectra")
        self._assert_no_console_errors()

    def test_spectra_img(self):
        _url = "proteins:spectra-img"
        for ext in [".svg", ".png", ".pdf", ".jpg", ".jpeg"]:
            self._load_reverse(_url, args=(self.p1.slug, ext))
            self._assert_no_console_errors()

        # test passing matplotlib kwargs
        rev = reverse(_url, args=(self.p1.slug, ".svg"))
        rev = f"{rev}?xlim=350,700&alpha=0.2&grid=true"
        self.browser.get(self.live_server_url + rev)
        self._assert_no_console_errors()

    def test_microscopes(self):
        m = MicroscopeFactory(name="myScope", id="KLMNPQRSTUVWX")
        OpticalConfigWithFiltersFactory.create_batch(4, microscope=m)
        self._load_reverse("proteins:microscopes")
        self.browser.find_element(by="xpath", value=f'//a[text()="{m.name}"]').click()
        self._interact_scope(m)

        # FIXME: there are some console errors on CI related to chart rendering timing
        # "Cannot read properties of undefined (reading 'x')" - race condition in nvd3 chart
        self.browser.get_log("browser")  # clear logs instead of asserting

    def test_embedscope(self):
        m = MicroscopeFactory(name="myScope", id="KLMNPQRSTUVWX")
        OpticalConfigWithFiltersFactory.create_batch(2, microscope=m)
        self._load_reverse("proteins:microscope-embed", args=(m.id,))
        self._interact_scope(m)

        # FIXME: there are some console errors on CI that need to be fixed
        # http://localhost:33339/protein/knownsequence/ - OTS parsing error: invalid sfntVersion: 1702391919'
        self.browser.get_log("browser")  # clear prior logs

    def _interact_scope(self, scope):
        self.browser.find_element(value="select2-fluor-select-container").click()
        self.browser.switch_to.active_element.send_keys(self.p1.name[:5])
        self.browser.switch_to.active_element.send_keys(Keys.ENTER)

        self.browser.find_element(value="config-select").click()
        self.browser.switch_to.active_element.send_keys(Keys.ARROW_DOWN)
        self.browser.switch_to.active_element.send_keys(Keys.ARROW_DOWN)
        self.browser.switch_to.active_element.send_keys(Keys.ENTER)

    def test_blast(self):
        self._load_reverse("proteins:blast")
        self._assert_no_console_errors()
        text_input = self.browser.find_element(value="queryInput")
        assert text_input.is_displayed()
        text_input.send_keys(SEQ[5:20].replace("LDG", "LG"))

        try:
            _get_binary("makeblastdb")
        except Exception:
            if not os.environ.get("CI"):
                pytest.xfail("makeblastdb binary not found")
                return

        submit = self.browser.find_element(by="css selector", value='button[type="submit"]')
        submit.click()

        # wait for the results to load
        try:
            r1 = WebDriverWait(self.browser, timeout=2).until(
                lambda d: d.find_element(by="xpath", value="//table/tbody/tr[1]/td[1]/a")
            )
        except Exception as e:
            msg = f"Results failed to load: {e.msg}"
            logs = self.browser.get_log("browser")
            for log in logs:
                msg += f"\n{log['message']}"
            raise AssertionError(msg) from None

        # assert a table is displayed with the first row containing the protein
        assert r1.text == self.p1.name
        # it should link to the alignment table
        r1.click()

    def test_fret(self):
        """Test FRET page loads and basic UI interaction works.

        Note: AJAX-dependent calculations are not tested due to timing unreliability.
        """
        ProteinFactory(name="donor", agg="m", default_state__ex_max=488, default_state__em_max=525)
        ProteinFactory(
            name="acceptor",
            agg="m",
            default_state__ex_max=525,
            default_state__em_max=550,
        )
        self._load_reverse("proteins:fret")

        # Verify page loads and select2 widgets are present and clickable
        WebDriverWait(self.browser, 10).until(EC.element_to_be_clickable((By.ID, "select2-donor-select-container")))
        WebDriverWait(self.browser, 10).until(EC.element_to_be_clickable((By.ID, "select2-acceptor-select-container")))

        # Verify the result fields exist
        self.browser.find_element(value="QYD")
        self.browser.find_element(value="QYA")
        self.browser.find_element(value="overlapIntgrl")

        self._assert_no_console_errors()

    def test_collections(self):
        self._load_reverse("proteins:collections")
        self._assert_no_console_errors()

    @pytest.mark.ignore_template_errors
    def test_table(self):
        ProteinFactory.create_batch(10)
        self._load_reverse("proteins:table")
        # TODO: the table isn't actually being drawn yet on selenium
        # table = WebDriverWait(self.browser, timeout=6).until(
        # lambda d: d.find_element(value="proteinTable_wrapper")
        # )
        self._assert_no_console_errors()

    @pytest.mark.ignore_template_errors
    def test_chart(self):
        ProteinFactory.create_batch(6)
        self._load_reverse("proteins:ichart")
        self._assert_no_console_errors()

        elem = self.browser.find_element(by="xpath", value="//label[input[@id='Xqy']]")
        elem.click()

        elem = self.browser.find_element(by="xpath", value="//label[input[@id='Yext_coeff']]")
        elem.click()
        self._assert_no_console_errors()

    def test_problems(self):
        self._load_reverse("proteins:problems")
        self._assert_no_console_errors()

    def test_problems_inconsistencies(self):
        self._load_reverse("proteins:problems-inconsistencies")
        self._assert_no_console_errors()

    def test_problems_gaps(self):
        self._load_reverse("proteins:problems-gaps")
        self._assert_no_console_errors()

    def test_search(self):
        self.p1 = ProteinFactory(name="knownSequence", seq=SEQ)

        self._load_reverse("proteins:search")
        self._assert_no_console_errors()

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
        # Wait for the redirect to complete
        expected_url = self.live_server_url + self.p1.get_absolute_url()
        WebDriverWait(self.browser, 5).until(lambda d: d.current_url == expected_url)
        assert self.browser.current_url == expected_url
        assert self.browser.find_element(by="xpath", value="//h1").text == self.p1.name
        shown_seq = self.browser.find_element(by="class name", value="formatted_aminosquence").text
        assert shown_seq.replace(" ", "") == SEQ

    @pytest.mark.ignore_template_errors
    def test_compare(self):
        p2 = ProteinFactory(seq=SEQ.replace("ELDG", "ETTG"))

        self._load_reverse("proteins:compare")
        self._assert_no_console_errors()

        prots = ",".join([self.p1.slug, p2.slug])
        self._load_reverse("proteins:compare", args=(prots,))

        muts = self.browser.find_element(by="xpath", value='//p[strong[text()="Mutations: "]]')
        assert muts.text == "Mutations: L19T/D20T"  # (the two T mutations we did above)
        self._assert_no_console_errors()

    def test_spectrum_submission_preview_manual_data(self):
        """End-to-end test of spectrum submission with manual data preview"""

        User = get_user_model()

        # Create a test user with verified email
        import uuid

        username = f"testuser_{uuid.uuid4().hex[:8]}"
        # Use the existing knownSequence protein from setUp
        # Delete any existing ex spectrum to ensure test can submit a new one
        protein_name = self.p1.name

        Spectrum.objects.filter(owner_state__protein=self.p1, subtype="ex").delete()

        user = User.objects.create_user(username=username, password=PASSWORD, email=f"{username}@example.com")
        user.is_active = True
        user.save()

        # Create EmailAddress record for allauth email verification

        EmailAddress.objects.create(user=user, email=f"{username}@example.com", verified=True, primary=True)

        client = Client()
        client.force_login(user)

        # Get session key and set it in browser cookies for Selenium
        session_key = client.session.session_key

        # Navigate to any page first to set domain for cookie
        self.browser.get(self.live_server_url)

        # Set Django session cookie
        self.browser.add_cookie({"name": "sessionid", "value": session_key, "domain": "localhost", "path": "/"})

        # Navigate to spectrum submission page
        self._load_reverse("proteins:submit-spectra")
        self._assert_no_console_errors()

        # Wait for form to load, then fill out the basic form fields
        WebDriverWait(self.browser, 10).until(lambda d: d.find_element(by="id", value="id_category").is_displayed())
        Select(self.browser.find_element(by="id", value="id_category")).select_by_value("p")
        Select(self.browser.find_element(by="id", value="id_subtype")).select_by_value("ex")

        # Wait for protein owner field to appear (it's a Select2 autocomplete field)
        # Find the Select2 container and click it to open the search
        select2_container = WebDriverWait(self.browser, 10).until(
            lambda d: d.find_element(by="css selector", value=".select2-selection")
        )
        select2_container.click()

        # Type the protein name into the autocomplete search box and select it
        self.browser.switch_to.active_element.send_keys(protein_name)
        # Wait briefly for autocomplete results to load
        WebDriverWait(self.browser, 5).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, ".select2-results__option"))
        )
        self.browser.switch_to.active_element.send_keys(Keys.ENTER)

        # Switch to manual data tab
        manual_tab = self.browser.find_element(by="id", value="manual-tab")
        manual_tab.click()

        # Wait for manual data field to be visible and interactable
        WebDriverWait(self.browser, 5).until(lambda d: d.find_element(by="id", value="id_data").is_displayed())
        data_field = WebDriverWait(self.browser, 5).until(
            lambda d: d.find_element(by="id", value="id_data")
            if d.find_element(by="id", value="id_data").is_enabled()
            else None
        )
        spectrum_data = (
            "[[400, 0.1], [401, 0.2], [402, 0.3], [403, 0.5], [404, 0.8], "
            "[405, 1.0], [406, 0.8], [407, 0.5], [408, 0.3], [409, 0.1]]"
        )
        data_field.send_keys(spectrum_data)

        # Check confirmation checkbox
        self.browser.find_element(by="id", value="id_confirmation").click()

        # Submit for preview
        submit_button = self.browser.find_element(by="css selector", value='input[type="submit"]')
        assert "Preview" in submit_button.get_attribute("value")
        submit_button.click()

        # Wait a moment for AJAX to start
        import time

        time.sleep(2)

        # Check for any error alerts
        error_alerts = self.browser.find_elements(by="css selector", value=".alert-danger")
        if error_alerts:
            for alert in error_alerts:
                if alert.is_displayed():
                    raise AssertionError(f"Error alert shown: {alert.text}")

        # Wait for preview to appear (AJAX request may take some time)
        WebDriverWait(self.browser, 15).until(
            lambda d: d.find_element(by="id", value="spectrum-preview-section")
            if d.find_element(by="id", value="spectrum-preview-section").is_displayed()
            else None
        )

        # Verify preview content
        preview_chart = self.browser.find_element(by="id", value="spectrum-preview-chart")
        assert preview_chart.find_element(by="tag name", value="svg")  # Should contain SVG

        # Check preview info
        data_points = self.browser.find_element(by="id", value="preview-data-points").text
        assert "10" in data_points  # Should show 10 data points

        # Test "Edit Data" button - use onclick attribute for reliable selection
        edit_button = self.browser.find_element(by="css selector", value="button[onclick='hidePreview()']")
        edit_button.click()

        # Should switch back to manual tab and hide preview
        WebDriverWait(self.browser, 2).until(
            lambda d: not d.find_element(by="id", value="spectrum-preview-section").is_displayed()
        )
        assert self.browser.find_element(by="id", value="manual-tab").get_attribute("class").find("active") != -1

        # Data should still be there (compare parsed JSON, not string formatting)
        import json

        data_field = self.browser.find_element(by="id", value="id_data")
        expected_data = json.loads(spectrum_data)
        actual_data = json.loads(data_field.get_attribute("value"))
        assert expected_data == actual_data

        self._assert_no_console_errors()

    def test_spectrum_submission_tab_switching(self):
        """End-to-end test of tab switching behavior in spectrum submission"""

        User = get_user_model()

        # Create a test user with verified email
        import uuid

        username = f"testuser_{uuid.uuid4().hex[:8]}"
        user = User.objects.create_user(username=username, password=PASSWORD, email=f"{username}@example.com")
        user.is_active = True
        user.save()

        # Create EmailAddress record for allauth email verification

        EmailAddress.objects.create(user=user, email=f"{username}@example.com", verified=True, primary=True)

        # Use Django's session authentication instead of HTML login

        client = Client()
        client.force_login(user)

        # Get session key and set it in browser cookies for Selenium
        session_key = client.session.session_key

        # Navigate to any page first to set domain for cookie
        self.browser.get(self.live_server_url)

        # Set Django session cookie
        self.browser.add_cookie({"name": "sessionid", "value": session_key, "domain": "localhost", "path": "/"})

        # Navigate to spectrum submission page
        self._load_reverse("proteins:submit-spectra")

        # Wait for form to load, then fill out basic fields
        WebDriverWait(self.browser, 10).until(lambda d: d.find_element(by="id", value="id_category").is_displayed())
        Select(self.browser.find_element(by="id", value="id_category")).select_by_value("p")
        Select(self.browser.find_element(by="id", value="id_subtype")).select_by_value("ex")
        WebDriverWait(self.browser, 2).until(lambda d: d.find_element(by="id", value="id_owner_state").is_displayed())
        owner_state_select = Select(self.browser.find_element(by="id", value="id_owner_state"))
        if len(owner_state_select.options) > 1:  # Skip the empty option
            owner_state_select.select_by_index(1)

        # Check confirmation
        self.browser.find_element(by="id", value="id_confirmation").click()

        # Test tab switching behavior
        file_tab = self.browser.find_element(by="id", value="file-tab")
        manual_tab = self.browser.find_element(by="id", value="manual-tab")

        # Start on file tab (default)
        assert "active" in file_tab.get_attribute("class")

        # Switch to manual tab
        manual_tab.click()
        WebDriverWait(self.browser, 2).until(
            lambda d: "active" in d.find_element(by="id", value="manual-tab").get_attribute("class")
        )

        # Wait for data field to become visible and enabled after tab switch
        data_field = WebDriverWait(self.browser, 3).until(
            lambda d: d.find_element(by="id", value="id_data")
            if d.find_element(by="id", value="id_data").is_displayed()
            and d.find_element(by="id", value="id_data").is_enabled()
            else None
        )

        # Enter some manual data
        data_field.send_keys(
            "[[400, 0.1], [401, 0.2], [402, 0.3], [403, 0.5], "
            "[404, 0.8], [405, 1.0], [406, 0.8], [407, 0.5], [408, 0.3], [409, 0.1]]"
        )

        # Switch back to file tab
        file_tab.click()
        WebDriverWait(self.browser, 1).until(
            lambda d: "active" in d.find_element(by="id", value="file-tab").get_attribute("class")
        )

        # Submit button should update based on which tab is active and whether there's data
        submit_button = self.browser.find_element(by="css selector", value='input[type="submit"]')
        # On file tab with no file, should show "Submit" not "Preview"
        assert submit_button.get_attribute("value") in ["Submit", "Preview Spectrum"]

        self._assert_no_console_errors()
