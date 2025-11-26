from django.test import RequestFactory, TestCase

from fpbase.users.views import UserRedirectView, UserUpdateView
from tests.test_users.factories import UserFactory


class BaseUserTestCase(TestCase):
    def setUp(self):
        self.user = UserFactory(username="testuser")
        self.factory = RequestFactory()


class TestUserRedirectView(BaseUserTestCase):
    def test_get_redirect_url(self):
        # Instantiate the view directly. Never do this outside a test!
        view = UserRedirectView()
        request = self.factory.get("/fake-url")
        request.user = self.user
        view.request = request


class TestUserUpdateView(BaseUserTestCase):
    def setUp(self):
        # call BaseUserTestCase.setUp()
        super().setUp()
        # Instantiate the view directly. Never do this outside a test!
        self.view = UserUpdateView()
        # Generate a fake request
        request = self.factory.get("/fake-url")
        # Attach the user to the request
        request.user = self.user
        # Attach the request to the view
        self.view.request = request

    def test_get_success_url(self):
        self.assertEqual(self.view.get_success_url(), f"/users/{self.user.username}/")

    def test_get_object(self):
        # Expect: self.user, as that is the request's user object
        self.assertEqual(self.view.get_object(), self.user)
