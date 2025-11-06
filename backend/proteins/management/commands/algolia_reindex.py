"""Management command to reindex Algolia."""

from django.core.management.base import BaseCommand
from proteins.tasks import reindex_all_proteins, configure_all_indices


class Command(BaseCommand):
    help = 'Reindex all models in Algolia'

    def add_arguments(self, parser):
        parser.add_argument(
            '--configure',
            action='store_true',
            help='Configure index settings before reindexing',
        )
        parser.add_argument(
            '--sync',
            action='store_true',
            help='Run synchronously (blocks until complete, useful for testing)',
        )

    def handle(self, *args, **options):
        if options['configure']:
            self.stdout.write('Configuring indices...')
            if options['sync']:
                # Run synchronously
                result = configure_all_indices()
                self.stdout.write(self.style.SUCCESS(f'✓ Indices configured: {result}'))
            else:
                # Queue as async task
                configure_all_indices.delay()
                self.stdout.write(self.style.SUCCESS('✓ Index configuration queued'))

        self.stdout.write('Reindexing proteins...')
        if options['sync']:
            # Run synchronously
            result = reindex_all_proteins()
            self.stdout.write(self.style.SUCCESS(f'✓ Reindexed {result["total_indexed"]} proteins'))
        else:
            # Queue as async task
            reindex_all_proteins.delay()
            self.stdout.write(self.style.SUCCESS('✓ Reindexing queued'))

        if not options['sync']:
            self.stdout.write('\nNote: Tasks are running asynchronously.')
            self.stdout.write('Check Celery logs for progress.')
