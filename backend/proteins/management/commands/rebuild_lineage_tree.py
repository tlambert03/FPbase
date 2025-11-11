from __future__ import annotations

from django.core.management.base import BaseCommand
from django.db import transaction

from proteins.models import Lineage


class Command(BaseCommand):
    help = "Rebuild the MPTT tree structure for Lineage model to fix any corruption"

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Show what would be done without making changes",
        )

    def handle(self, **options):
        dry_run = options["dry_run"]

        self.stdout.write("Analyzing current Lineage tree structure...")

        # Get statistics before rebuild
        total_lineages = Lineage.objects.count()
        tree_ids = Lineage.objects.values_list("tree_id", flat=True).distinct()
        num_trees = len(set(tree_ids))

        self.stdout.write(f"  Total lineage nodes: {total_lineages}")
        self.stdout.write(f"  Number of distinct trees: {num_trees}")

        # Check for potential corruption: nodes with parents in different trees
        corrupt_nodes = []
        for lineage in Lineage.objects.select_related("parent").exclude(parent=None):
            if lineage.parent and lineage.tree_id != lineage.parent.tree_id:
                corrupt_nodes.append(
                    f"    - {lineage.protein.name} (ID {lineage.id}, "
                    f"tree_id={lineage.tree_id}) has parent "
                    f"{lineage.parent.protein.name} (ID {lineage.parent.id}, "
                    f"tree_id={lineage.parent.tree_id})"
                )

        if corrupt_nodes:
            self.stdout.write(self.style.WARNING(f"\nFound {len(corrupt_nodes)} nodes with tree_id inconsistencies:"))
            for node in corrupt_nodes[:10]:  # Show first 10
                self.stdout.write(self.style.WARNING(node))
            if len(corrupt_nodes) > 10:
                self.stdout.write(self.style.WARNING(f"    ... and {len(corrupt_nodes) - 10} more"))
        else:
            self.stdout.write(self.style.SUCCESS("\n  No tree_id inconsistencies detected."))

        if dry_run:
            self.stdout.write("\n" + self.style.WARNING("DRY RUN: No changes will be made."))
            return

        self.stdout.write("\nRebuilding MPTT tree structure...")

        try:
            with transaction.atomic():
                Lineage.objects.rebuild()

            # Get statistics after rebuild
            new_tree_ids = Lineage.objects.values_list("tree_id", flat=True).distinct()
            new_num_trees = len(set(new_tree_ids))

            self.stdout.write(self.style.SUCCESS("\n✓ MPTT tree rebuild completed successfully!"))
            self.stdout.write(f"\n  Trees before: {num_trees}")
            self.stdout.write(f"  Trees after:  {new_num_trees}")

            if corrupt_nodes:
                self.stdout.write(self.style.SUCCESS(f"\n  Fixed {len(corrupt_nodes)} tree_id inconsistencies"))

        except Exception as e:
            self.stdout.write(self.style.ERROR(f"\n✗ Rebuild failed: {e}"))
            raise
