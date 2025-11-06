"""
Django management command for analyzing memory usage.
Run with: uv run backend/manage.py memory_analysis
"""
import gc
import sys

from django.core.management.base import BaseCommand
from pympler import asizeof, muppy, summary


class Command(BaseCommand):
    help = "Analyze memory usage of Django application"

    def add_arguments(self, parser):
        parser.add_argument(
            "--top",
            type=int,
            default=20,
            help="Number of top objects to show (default: 20)",
        )
        parser.add_argument(
            "--type",
            type=str,
            help="Filter by object type (e.g., 'dict', 'list')",
        )

    def handle(self, *args, **options):
        self.stdout.write(self.style.SUCCESS("Starting memory analysis..."))

        # Force garbage collection
        gc.collect()

        # Get all objects
        all_objects = muppy.get_objects()

        # Get summary
        sum1 = summary.summarize(all_objects)

        # Print summary
        self.stdout.write("\n" + "=" * 80)
        self.stdout.write("MEMORY USAGE SUMMARY")
        self.stdout.write("=" * 80 + "\n")

        summary.print_(sum1, limit=options["top"])

        # Filter by type if specified
        if options["type"]:
            self.stdout.write(f"\n\nFiltering for type: {options['type']}")
            filtered = [obj for obj in all_objects if type(obj).__name__ == options["type"]]
            self.stdout.write(f"Found {len(filtered)} objects of type {options['type']}")
            self.stdout.write(f"Total size: {asizeof.asizeof(filtered) / 1024 / 1024:.2f} MB")

        # Django-specific analysis
        self.stdout.write("\n" + "=" * 80)
        self.stdout.write("DJANGO-SPECIFIC ANALYSIS")
        self.stdout.write("=" * 80 + "\n")

        # Check for common Django memory issues
        from django.db import connection

        self.stdout.write(f"Database queries: {len(connection.queries)}")
        self.stdout.write(
            f"Query memory: {asizeof.asizeof(connection.queries) / 1024:.2f} KB"
        )

        # Check Django settings
        from django.conf import settings

        self.stdout.write(f"Settings size: {asizeof.asizeof(settings._wrapped) / 1024:.2f} KB")

        self.stdout.write("\n" + self.style.SUCCESS("Memory analysis complete!"))
