from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model
from django.test import TestCase

from proteins.factories import ProteinFactory

from .models import Favorite

User = get_user_model()


@pytest.mark.django_db
class TestFavoriteRaceCondition:
    """Test that the favorite creation handles race conditions correctly."""

    def test_get_or_create_handles_concurrent_creation(self):
        """Test that get_or_create prevents duplicate favorites."""
        user = User.objects.create_user(username="testuser", password="testpass")
        protein = ProteinFactory()

        # Try to create the same favorite twice using get_or_create
        fav1, created1 = Favorite.objects.get_or_create(user, protein.id, "proteins.Protein")
        fav2, created2 = Favorite.objects.get_or_create(user, protein.id, "proteins.Protein")

        # First should be created, second should be retrieved
        assert created1 is True
        assert created2 is False
        assert fav1.id == fav2.id

        # Should only have one favorite in database
        assert Favorite.objects.filter(user=user).count() == 1

    def test_multiple_get_or_create_calls_return_same_favorite(self):
        """Test that multiple get_or_create calls return the same favorite."""
        user = User.objects.create_user(username="testuser", password="testpass")
        protein = ProteinFactory()

        # Call get_or_create multiple times
        results = []
        for _ in range(5):
            fav, created = Favorite.objects.get_or_create(user, protein.id, "proteins.Protein")
            results.append((fav.id, created))

        # All should return the same favorite ID
        fav_ids = [fav_id for fav_id, _ in results]
        assert len(set(fav_ids)) == 1, "All calls should return the same favorite ID"

        # Only first should be marked as created
        created_count = sum(1 for _, created in results if created)
        assert created_count == 1, "Only first call should have created the favorite"

        # Verify only one favorite exists in database
        assert Favorite.objects.filter(user=user).count() == 1


class TestFavoriteManager(TestCase):
    """Test FavoriteManager methods."""

    def setUp(self):
        self.user = User.objects.create_user(username="testuser", password="testpass")
        self.protein = ProteinFactory()

    def test_get_favorite_returns_none_when_not_exists(self):
        """Test that get_favorite returns None when favorite doesn't exist."""
        fav = Favorite.objects.get_favorite(self.user, self.protein.id, "proteins.Protein")
        assert fav is None

    def test_get_favorite_returns_favorite_when_exists(self):
        """Test that get_favorite returns favorite when it exists."""
        created_fav = Favorite.objects.create(self.user, self.protein.id, "proteins.Protein")
        retrieved_fav = Favorite.objects.get_favorite(self.user, self.protein.id, "proteins.Protein")

        assert retrieved_fav is not None
        assert retrieved_fav.id == created_fav.id

    def test_for_user_returns_user_favorites(self):
        """Test that for_user returns all favorites for a user."""
        protein2 = ProteinFactory()

        Favorite.objects.create(self.user, self.protein.id, "proteins.Protein")
        Favorite.objects.create(self.user, protein2.id, "proteins.Protein")

        favorites = Favorite.objects.for_user(self.user)
        assert favorites.count() == 2

    def test_for_object_returns_object_favorites(self):
        """Test that for_object returns all favorites for an object."""
        user2 = User.objects.create_user(username="testuser2", password="testpass")

        Favorite.objects.create(self.user, self.protein.id, "proteins.Protein")
        Favorite.objects.create(user2, self.protein.id, "proteins.Protein")

        favorites = Favorite.objects.for_object(self.protein.id, "proteins.Protein")
        assert favorites.count() == 2
