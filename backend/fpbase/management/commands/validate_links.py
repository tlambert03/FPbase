"""Management command to validate outgoing links."""

from django.core.management.base import BaseCommand

from fpbase.tasks import validate_outgoing_links


class Command(BaseCommand):
    help = "Validate all outgoing links on the site"

    def add_arguments(self, parser):
        parser.add_argument(
            "--sources",
            nargs="+",
            choices=["database", "templates", "static"],
            help="Link sources to check (default: all)",
        )
        parser.add_argument(
            "--async",
            action="store_true",
            dest="run_async",
            help="Run as a background Celery task",
        )
        parser.add_argument(
            "--fix",
            action="store_true",
            help="Attempt to fix broken links (not yet implemented)",
        )

    def handle(self, *args, **options):
        sources = options.get("sources")
        run_async = options.get("run_async", False)
        fix_broken = options.get("fix", False)

        self.stdout.write("Starting link validation...")

        if run_async:
            # Run as Celery task
            task = validate_outgoing_links.delay(sources=sources, fix_broken=fix_broken)
            self.stdout.write(
                self.style.SUCCESS(
                    f"Link validation task queued with ID: {task.id}\n"
                    f"Check status with: celery -A fpbase.celery inspect active"
                )
            )
        else:
            # Run synchronously
            results = validate_outgoing_links(sources=sources, fix_broken=fix_broken)

            # Display results
            self.stdout.write("\n" + "=" * 70)
            self.stdout.write(self.style.SUCCESS("Link Validation Results"))
            self.stdout.write("=" * 70)
            self.stdout.write(f"Total links checked: {results['total_links']}")
            self.stdout.write(
                self.style.SUCCESS(f"Valid links: {results['valid_links']}")
            )
            self.stdout.write(
                self.style.WARNING(f"Skipped links: {results['skipped_links']}")
            )
            self.stdout.write(
                self.style.ERROR(f"Broken links: {len(results['broken_links'])}")
            )

            if results["broken_links"]:
                self.stdout.write("\n" + self.style.ERROR("Broken Links:"))
                for broken in results["broken_links"]:
                    self.stdout.write(f"  - URL: {broken.get('url', 'N/A')}")
                    self.stdout.write(f"    Source: {broken.get('source', 'N/A')}")
                    self.stdout.write(f"    Error: {broken.get('error', 'N/A')}")
                    if "context" in broken:
                        self.stdout.write(f"    Context: {broken['context']}")
                    self.stdout.write("")

            self.stdout.write("=" * 70 + "\n")
