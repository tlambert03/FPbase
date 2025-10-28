from django.test import TestCase

from .factories import UserFactory


class TestUser(TestCase):
    def setUp(self):
        self.user = UserFactory()

    def test__str__(self):
        self.assertEqual(self.user.__str__(), self.user.username)

    def test_get_absolute_url(self):
        self.assertEqual(self.user.get_absolute_url(), f"/users/{self.user.username}/")
