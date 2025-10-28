from django.test import TestCase
from django.urls import resolve, reverse

from .factories import UserFactory


class TestUserURLs(TestCase):
    """Test URL patterns for users app."""

    def setUp(self):
        self.uname = "testuser"
        self.user = UserFactory(username=self.uname)

    def test_list_reverse(self):
        """users:list should reverse to /users/."""
        self.assertEqual(reverse("users:list"), "/users/")

    def test_list_resolve(self):
        """/users/ should resolve to users:list."""
        self.assertEqual(resolve("/users/").view_name, "users:list")

    def test_redirect_reverse(self):
        """users:redirect should reverse to /users/~redirect/."""
        self.assertEqual(reverse("users:redirect"), "/users/~redirect/")

    def test_redirect_resolve(self):
        """/users/~redirect/ should resolve to users:redirect."""
        self.assertEqual(resolve("/users/~redirect/").view_name, "users:redirect")

    def test_detail_reverse(self):
        """users:detail should reverse to /users/testuser/."""
        self.assertEqual(reverse("users:detail", kwargs={"username": f"{self.uname}"}), f"/users/{self.uname}/")

    def test_detail_resolve(self):
        """/users/testuser/ should resolve to users:detail."""
        self.assertEqual(resolve(f"/users/{self.uname}/").view_name, "users:detail")

    def test_update_reverse(self):
        """users:update should reverse to /users/~update/."""
        self.assertEqual(reverse("users:update"), "/users/~update/")

    def test_update_resolve(self):
        """/users/~update/ should resolve to users:update."""
        self.assertEqual(resolve("/users/~update/").view_name, "users:update")
