from __future__ import annotations

import pytest

from proteins.models import Lineage, Protein


@pytest.mark.django_db
class TestLineageMPTT:
    """Test django-mptt functionality in Lineage model."""

    def test_root_lineage_creation(self):
        """Test creating a root lineage node sets MPTT fields correctly."""
        protein = Protein.objects.create(name="TestProtein", slug="test-protein")
        lineage = Lineage.objects.create(protein=protein)

        # MPTT fields should be populated
        assert lineage.lft is not None
        assert lineage.rght is not None
        assert lineage.tree_id is not None
        assert lineage.level == 0
        assert lineage.get_root() == lineage

    def test_child_lineage_tree_structure(self):
        """Test parent-child relationships maintain proper tree structure."""
        # Create root
        root_protein = Protein.objects.create(name="RootProtein", slug="root")
        root_lineage = Lineage.objects.create(protein=root_protein)

        # Create child
        child_protein = Protein.objects.create(name="ChildProtein", slug="child")
        child_lineage = Lineage.objects.create(protein=child_protein, parent=root_lineage)

        # Verify tree structure
        assert child_lineage.parent == root_lineage
        assert child_lineage.get_root() == root_lineage
        assert child_lineage.level == 1
        assert child_lineage in root_lineage.get_children()
        assert child_lineage in root_lineage.get_descendants()

        # Create grandchild
        grandchild_protein = Protein.objects.create(name="GrandchildProtein", slug="grandchild")
        grandchild_lineage = Lineage.objects.create(
            protein=grandchild_protein, parent=child_lineage
        )

        # Verify multi-level tree
        assert grandchild_lineage.get_root() == root_lineage
        assert grandchild_lineage.level == 2
        assert list(grandchild_lineage.get_ancestors()) == [
            root_lineage,
            child_lineage,
        ]
        assert grandchild_lineage in root_lineage.get_descendants()
        assert root_lineage.get_descendant_count() == 2

    def test_multiple_trees(self):
        """Test multiple root nodes create separate trees."""
        # Create two separate trees
        protein1 = Protein.objects.create(name="Tree1Root", slug="tree1")
        lineage1 = Lineage.objects.create(protein=protein1)

        protein2 = Protein.objects.create(name="Tree2Root", slug="tree2")
        lineage2 = Lineage.objects.create(protein=protein2)

        # Verify they're in different trees
        assert lineage1.tree_id != lineage2.tree_id
        assert lineage2 not in lineage1.get_descendants()
        assert lineage1 not in lineage2.get_descendants()
