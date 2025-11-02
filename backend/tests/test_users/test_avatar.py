"""Tests for django-avatar functionality."""

import shutil
from io import BytesIO

import pytest
from avatar.models import Avatar
from avatar.templatetags.avatar_tags import avatar_url
from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse
from PIL import Image

from .factories import UserFactory

User = get_user_model()


def create_test_image():
    """Create a simple test image for avatar uploads."""
    image = Image.new("RGB", (100, 100), color="red")
    img_io = BytesIO()
    image.save(img_io, format="PNG")
    img_io.seek(0)
    return SimpleUploadedFile("test_avatar.png", img_io.read(), content_type="image/png")


@pytest.fixture
def temp_media_root(tmp_path):
    """Create a temporary media root for tests."""
    media_root = tmp_path / "test_media"
    media_root.mkdir()
    yield media_root
    # Cleanup happens automatically with tmp_path


class AvatarTestCase(TestCase):
    """Base test case with authenticated user and temporary media storage."""

    @pytest.fixture(autouse=True)
    def setup_temp_media(self, tmp_path):
        """Set up temporary media directory for each test."""
        self.temp_media = tmp_path / "test_media"
        self.temp_media.mkdir()
        # Store original MEDIA_ROOT to restore after test
        from django.conf import settings

        self._original_media_root = settings.MEDIA_ROOT
        settings.MEDIA_ROOT = str(self.temp_media)

    def setUp(self):
        self.user = UserFactory()
        self.client.force_login(self.user)

    def tearDown(self):
        """Clean up test media files."""
        # Restore original MEDIA_ROOT
        from django.conf import settings

        settings.MEDIA_ROOT = self._original_media_root
        # Clean up any created files
        if hasattr(self, "temp_media") and self.temp_media.exists():
            shutil.rmtree(self.temp_media, ignore_errors=True)


class TestAvatarUpload(AvatarTestCase):
    """Test avatar upload functionality."""

    def test_avatar_upload_view_accessible(self):
        """Test that avatar add page is accessible."""
        url = reverse("avatar:add")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)

    def test_avatar_upload_requires_authentication(self):
        """Test that avatar upload requires login."""
        self.client.logout()
        url = reverse("avatar:add")
        response = self.client.get(url)
        # Should redirect to login
        self.assertEqual(response.status_code, 302)
        self.assertIn("/accounts/login/", response.url)

    @override_settings(AVATAR_AUTO_GENERATE_SIZES=(80, 40))
    def test_avatar_upload_creates_avatar(self):
        """Test that uploading an image creates an avatar."""
        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 0)

        url = reverse("avatar:add")
        image = create_test_image()
        response = self.client.post(url, {"avatar": image}, follow=True)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 1)

        avatar = Avatar.objects.get(user=self.user)
        self.assertTrue(avatar.primary)
        self.assertTrue(avatar.avatar.name.endswith(".png"))

    def test_avatar_upload_invalid_file(self):
        """Test that uploading non-image file fails gracefully."""
        url = reverse("avatar:add")
        invalid_file = SimpleUploadedFile("test.txt", b"not an image", content_type="text/plain")
        response = self.client.post(url, {"avatar": invalid_file})

        # Should show form with errors
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.context["upload_avatar_form"].errors)
        self.assertIn("avatar", response.context["upload_avatar_form"].errors)


class TestAvatarChange(AvatarTestCase):
    """Test avatar change (selection) functionality."""

    def setUp(self):
        super().setUp()
        # Create two avatars for the user
        self.avatar1 = Avatar.objects.create(
            user=self.user, primary=True, avatar=SimpleUploadedFile("avatar1.png", b"image1", content_type="image/png")
        )
        self.avatar2 = Avatar.objects.create(
            user=self.user,
            primary=False,
            avatar=SimpleUploadedFile("avatar2.png", b"image2", content_type="image/png"),
        )

    def test_avatar_change_view_accessible(self):
        """Test that avatar change page is accessible."""
        url = reverse("avatar:change")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Your current avatar:")

    def test_avatar_change_switches_primary(self):
        """Test that changing avatar switches the primary flag."""
        url = reverse("avatar:change")
        response = self.client.post(url, {"choice": str(self.avatar2.id)}, follow=True)

        self.assertEqual(response.status_code, 200)

        # Refresh from database
        self.avatar1.refresh_from_db()
        self.avatar2.refresh_from_db()

        self.assertFalse(self.avatar1.primary)
        self.assertTrue(self.avatar2.primary)

    def test_avatar_change_requires_authentication(self):
        """Test that avatar change requires login."""
        self.client.logout()
        url = reverse("avatar:change")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 302)
        self.assertIn("/accounts/login/", response.url)


class TestAvatarDelete(AvatarTestCase):
    """Test avatar deletion functionality."""

    def setUp(self):
        super().setUp()
        self.avatar = Avatar.objects.create(
            user=self.user, primary=True, avatar=SimpleUploadedFile("avatar.png", b"image", content_type="image/png")
        )

    def test_avatar_delete_view_accessible(self):
        """Test that avatar delete page is accessible."""
        url = reverse("avatar:delete")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Please select the avatars that you would like to delete")

    def test_avatar_delete_removes_avatar(self):
        """Test that deleting avatar removes it from database."""
        url = reverse("avatar:delete")
        response = self.client.post(url, {"choices": [str(self.avatar.id)]}, follow=True)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 0)

    def test_avatar_delete_requires_authentication(self):
        """Test that avatar delete requires login."""
        self.client.logout()
        url = reverse("avatar:delete")
        response = self.client.get(url)
        self.assertEqual(response.status_code, 302)
        self.assertIn("/accounts/login/", response.url)


class TestAvatarTemplateTag(AvatarTestCase):
    """Test avatar_url template tag."""

    def setUp(self):
        super().setUp()
        self.avatar = Avatar.objects.create(
            user=self.user, primary=True, avatar=SimpleUploadedFile("avatar.png", b"image", content_type="image/png")
        )

    def test_avatar_url_returns_url(self):
        """Test that avatar_url returns a valid URL."""
        url = avatar_url(self.user)
        self.assertIsNotNone(url)
        self.assertIn("/media/avatars/", url)

    def test_avatar_url_with_size(self):
        """Test that avatar_url accepts size parameter."""
        url = avatar_url(self.user, 80)
        self.assertIsNotNone(url)

    def test_avatar_url_without_avatar(self):
        """Test avatar_url for user without avatar returns default."""
        user_without_avatar = UserFactory()
        url = avatar_url(user_without_avatar)
        # Should return a default/gravatar URL
        self.assertIsNotNone(url)


class TestAvatarIntegration(AvatarTestCase):
    """Integration tests for avatar functionality."""

    def test_user_detail_page_shows_avatar_link(self):
        """Test that user detail page shows change avatar link."""
        url = reverse("users:detail", kwargs={"username": self.user.username})
        response = self.client.get(url, headers={"host": "testserver"})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Change your avatar")
        self.assertContains(response, reverse("avatar:change"))

    def test_avatar_workflow_upload_change_delete(self):
        """Test complete avatar workflow: upload, change, delete."""
        # 1. Upload first avatar
        url = reverse("avatar:add")
        image1 = create_test_image()
        response = self.client.post(url, {"avatar": image1}, follow=True)
        self.assertEqual(response.status_code, 200)
        avatar1_id = Avatar.objects.get(user=self.user).id

        # 2. Upload second avatar
        image2 = create_test_image()
        response = self.client.post(url, {"avatar": image2}, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 2)

        # Find the second avatar
        avatar2 = Avatar.objects.filter(user=self.user).exclude(id=avatar1_id).first()

        # 3. Change primary avatar
        url = reverse("avatar:change")
        response = self.client.post(url, {"choice": str(avatar2.id)}, follow=True)
        self.assertEqual(response.status_code, 200)
        avatar2.refresh_from_db()
        self.assertTrue(avatar2.primary)

        # 4. Delete avatars
        url = reverse("avatar:delete")
        response = self.client.post(url, {"choices": [str(avatar1_id), str(avatar2.id)]}, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 0)

    def test_multiple_users_can_have_avatars(self):
        """Test that multiple users can each have their own avatars."""
        user2 = UserFactory()

        # Create avatar for user1
        avatar1 = Avatar.objects.create(
            user=self.user, primary=True, avatar=SimpleUploadedFile("u1.png", b"img1", content_type="image/png")
        )

        # Create avatar for user2
        avatar2 = Avatar.objects.create(
            user=user2, primary=True, avatar=SimpleUploadedFile("u2.png", b"img2", content_type="image/png")
        )

        self.assertEqual(Avatar.objects.filter(user=self.user).count(), 1)
        self.assertEqual(Avatar.objects.filter(user=user2).count(), 1)
        self.assertNotEqual(avatar1.avatar.name, avatar2.avatar.name)
