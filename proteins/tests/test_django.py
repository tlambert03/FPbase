# from fpbase.users.models import User
# from django.contrib.staticfiles.testing import StaticLiveServerTestCase
# from allauth.account.models import EmailAddress
# from selenium import webdriver
# from selenium.webdriver.common.by import By
# from selenium.webdriver.support.ui import WebDriverWait
# from selenium.webdriver.support import expected_conditions as EC


# class UserRegistrationSeleniumTestCase(StaticLiveServerTestCase):

#     def setUp(self):
#         self.browser = webdriver.Chrome()
#         self.browser.get(self.live_server_url)

#     def test_user_registration(self):
#         self.browser.find_element_by_id("log-in-link").click()
#         self.browser.find_element_by_id("signup-link").click()

#         username = "newuser"
#         self.browser.find_element_by_id("id_username").send_keys(username)
#         self.browser.find_element_by_id("id_email").send_keys("newuser@email.com")
#         self.browser.find_element_by_id("id_password1").send_keys("Psiph5sK")
#         self.browser.find_element_by_id("id_password2").send_keys("Psiph5sK")
#         self.browser.find_element_by_id("id_password2").submit()

#         #self.assertEqual(username, self.browser.find_element_by_id("username-text").text)


# class UserLoginSeleniumTestCase(StaticLiveServerTestCase):

#     def setUp(self):
#         self.browser = webdriver.Chrome()
#         self.browser.get(self.live_server_url)
#         email = "todo@todoapp.com"
#         self.user = User.objects.create_user(username="newuser", password="NiGiw3Ch", email=email)
#         emails = EmailAddress.objects.filter(email=email, verified=False)
#         for email in emails:
#             email.verified = True
#             email.save()

#     def tearDown(self):
#         self.browser.quit()

#     def test_user_login(self):
#         self.browser.find_element_by_id("log-in-link").click()
#         self.browser.find_element_by_id("id_login").send_keys("newuser")
#         self.browser.find_element_by_id("id_password").send_keys("NiGiw3Ch")
#         self.browser.find_element_by_id("id_password").submit()
#         try:
#             element = WebDriverWait(self.browser, 10).until(
#                 EC.presence_of_element_located((By.CSS_SELECTOR, "div.alert-success"))
#             )
#             self.assertEqual('Successfully signed in as {}.'.format(self.user.username),
#                 element.text)
#         except Exception:
#             raise
