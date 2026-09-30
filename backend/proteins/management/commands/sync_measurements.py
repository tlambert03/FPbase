"""Record fluorescence values that live only on a state as measurements."""

from __future__ import annotations

from django.core.management.base import BaseCommand
from django.db import transaction

from proteins.models import DyeState, State


class Command(BaseCommand):
    help = (
        "Find states whose fluorescence values differ from what their measurements say "
        "(e.g. edited before edits were written through to measurements) and record those "
        "values as a measurement for the owner's primary reference. Idempotent."
    )

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true", help="only report")

    @transaction.atomic
    def handle(self, *args, dry_run: bool = False, **options):
        count = 0
        for model, owner in ((State, "protein"), (DyeState, "dye")):
            for state in model.objects.select_related(owner).iterator():
                composite, _ = state._composite()
                fields = [
                    f for f in state._written_through_fields() if getattr(state, f) != composite[f]
                ]
                if not fields:
                    continue
                count += 1
                changes = ", ".join(f"{f}: {composite[f]} -> {getattr(state, f)}" for f in fields)
                self.stdout.write(f"{state.slug}: {changes}")
                if not dry_run:
                    state.write_through(fields)
        verb = "would be" if dry_run else "were"
        self.stdout.write(f"{count} state(s) {verb} synced")
