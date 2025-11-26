from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING

import httpx
from django.core.management.base import BaseCommand
from django.db import transaction

from proteins.models import Protein, SnapGenePlasmid

if TYPE_CHECKING:
    from collections.abc import Sequence


class Command(BaseCommand):
    help = "Fetch and sync SnapGene plasmid data from their website"
    """
    This command fetches plasmid data from SnapGene's fluorescent protein plasmids page
    and syncs it to the database, then matches plasmids to proteins.

    Usage:
        python manage.py sync_snapgene_plasmids         # Run sync
        python manage.py sync_snapgene_plasmids --dry-run  # Preview without changes

    To set up as a periodic task (recommended to run monthly):
    - Add to Celery beat schedule, or
    - Set up as a cron job, or
    - Run manually when needed
    """

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Show what would be done without making changes",
        )

    def handle(self, **options):
        dry_run = options["dry_run"]

        if not (plasmid_data := self._fetch_plasmid_data()):
            return

        self.stdout.write(f"Found {len(plasmid_data)} plasmids")

        if dry_run:
            self.stdout.write("\n" + self.style.WARNING("DRY RUN: No changes will be made."))
            self._preview_plasmids(plasmid_data[:20])
            return

        # Sync plasmids to database
        created_count, updated_count, unchanged_count, removed_count = self._sync_plasmids(
            plasmid_data
        )
        self.stdout.write(
            self.style.SUCCESS(
                f"\n✓ Synced plasmids: {created_count} created, {updated_count} updated, "
                f"{unchanged_count} unchanged, {removed_count} removed"
            )
        )

        # Match plasmids to proteins
        matched_count = self._match_plasmids_to_proteins()
        self.stdout.write(self.style.SUCCESS(f"✓ Matched {matched_count} plasmids to proteins"))

    def _fetch_plasmid_data(self) -> list[dict]:
        self.stdout.write("Fetching SnapGene plasmid data...")
        try:
            # Fetch the page containing the plasmid data
            response = httpx.get(
                "https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids",
                timeout=30.0,
            )
            response.raise_for_status()
        except httpx.HTTPError as e:
            self.stdout.write(self.style.ERROR(f"✗ Failed to fetch data: {e}"))
            return []

        data = self._extract_plasmid_data(response.text)
        if not data:
            self.stdout.write(self.style.ERROR("✗ No plasmid data found in response"))
        return data

    def _extract_plasmid_data(self, content: str) -> list[dict]:
        """Extract plasmid data from the HTML page."""
        # Look for the JSON data embedded in the page
        # The data is in a variable like: var plasmidData = {...}
        match = re.search(r"var\s+plasmidData\s*=\s*({.*?});", content, re.DOTALL)
        if not match:
            # Try alternative pattern
            match = re.search(r'"sequences"\s*:\s*\[(.*?)\]', content, re.DOTALL)
            if not match:
                return []

        try:
            data = json.loads(
                match.group(1) if match.group(0).startswith("var") else f"[{match.group(1)}]"
            )
            # Handle different possible structures
            if isinstance(data, dict) and "sequences" in data:
                return data["sequences"]
            if isinstance(data, list):
                return data
        except json.JSONDecodeError as e:
            self.stdout.write(self.style.ERROR(f"✗ Failed to parse plasmid data: {e}"))
            pass

        return []

    def _preview_plasmids(self, plasmids: Sequence[dict]) -> None:
        """Show preview of plasmids that would be synced."""
        self.stdout.write("\nPreview of plasmids (first 10):")
        for p in plasmids:
            self.stdout.write(
                f"  - {p.get('plasmidName', 'Unknown')} ({p.get('plasmidID', 'Unknown')})"
            )

    def _sync_plasmids(self, plasmid_data: list[dict]) -> tuple[int, int, int, int]:
        """Sync plasmid data using get() then create-only logic.

        Tries to `get()` existing rows; if found counts as unchanged (no update).
        If missing, uses `update_or_create` (to satisfy requested pattern) to create.
        Removes plasmids that are no longer in the fetched data.
        Returns (created_count, updated_count, unchanged_count, removed_count).
        Updated will normally be 0 unless a race occurred between get and update_or_create.
        """
        created_count = 0
        updated_count = 0  # retained for signature; normally 0
        unchanged_count = 0
        removed_count = 0

        with transaction.atomic():
            # Track which plasmid IDs we've seen from SnapGene
            fetched_plasmid_ids = set()

            for data in plasmid_data:
                plasmid_id = data.get("plasmidID")
                if not plasmid_id:
                    self.stdout.write(
                        self.style.WARNING(
                            "✗ no plasmidID found... has the page structure changed?"
                        )
                    )
                    continue

                fetched_plasmid_ids.add(plasmid_id)

                defaults = {
                    "name": data.get("plasmidName", ""),
                    "description": data.get("plasmidDesc", ""),
                    "author": data.get("plasmidAuthor", ""),
                    "size": data.get("plasmidSize"),
                    "topology": data.get("plasmidTopology", ""),
                }

                try:
                    # Unchanged only if ALL fields already match current defaults.
                    SnapGenePlasmid.objects.get(plasmid_id=plasmid_id, **defaults)
                    unchanged_count += 1
                except SnapGenePlasmid.DoesNotExist:
                    _, created = SnapGenePlasmid.objects.update_or_create(
                        plasmid_id=plasmid_id, defaults=defaults
                    )
                    if created:
                        created_count += 1
                    else:
                        updated_count += 1

            # Remove plasmids that are no longer on SnapGene's website
            removed_count, _ = SnapGenePlasmid.objects.exclude(
                plasmid_id__in=fetched_plasmid_ids
            ).delete()

        return created_count, updated_count, unchanged_count, removed_count

    def _match_plasmids_to_proteins(self) -> int:
        """Match plasmids to proteins based on naming patterns."""
        matched_count = 0

        with transaction.atomic():
            # Clear existing relationships
            for protein in Protein.objects.all():
                protein.snapgene_plasmids.clear()

            # Match plasmids to proteins
            for protein in Protein.objects.all():
                protein_name = protein.name

                # Pattern to match exact name or plasmid vectors (e.g., pEGFP, pEGFP-N1, pEGFP-C2)
                is_plasmid_vector = re.compile(
                    rf"^p{re.escape(protein_name)}(-[CN]?\d+)?$", re.IGNORECASE
                )
                matching_plasmids = [
                    plasmid
                    for plasmid in SnapGenePlasmid.objects.all()
                    if plasmid.name == protein_name or is_plasmid_vector.match(plasmid.name)
                ]

                if matching_plasmids:
                    protein.snapgene_plasmids.set(matching_plasmids)
                    matched_count += len(matching_plasmids)

        return matched_count
