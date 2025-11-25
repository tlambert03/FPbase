"""Major data migration to new Fluorophore + Measurement schema.

NON REVERSIBLE.

1. proteins.State becomes MTI child of new proteins.Fluorophore
2. proteins.Dye becomes container for new proteins.DyeState (MTI child of Fluorophore)
3. proteins.FluorescenceMeasurement created for each old State and DyeState
4. proteins.Spectrum ownership updated to point to new Fluorophore records
5. proteins.OcFluorEff updated to point to new Fluorophore records

"""

from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING, Any

import django.db.models.deletion
import django.utils.timezone
import model_utils.fields
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import migrations, models

if TYPE_CHECKING:
    from django.apps.registry import Apps
    from django.db.backends.base.schema import BaseDatabaseSchemaEditor
    from django.db.backends.utils import CursorWrapper

logger = logging.getLogger(__name__)


def _dictfetchall(cursor: CursorWrapper) -> list[dict[str, Any]]:
    """Return all rows from a cursor as a dict. Assume the column names are unique."""
    if not cursor.description:
        return []
    columns = [col.name for col in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


MEASURABLE_FIELDS = {
    "ex_max",
    "em_max",
    "ext_coeff",
    "qy",
    "brightness",
    "lifetime",
    "pka",
    "twop_ex_max",
    "twop_peakGM",
    "twop_qy",
    "is_dark",
    "emhex",
    "exhex",
}


def migrate_state_data(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Migrate State data from old schema to new Fluorophore + State MTI structure."""
    # Get models from migration state
    State = apps.get_model("proteins", "State")
    Protein = apps.get_model("proteins", "Protein")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    # Access old State data directly via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy,
                   maturation, protein_id, created_by_id, updated_by_id, emhex, exhex
            FROM proteins_state_old
        """)
        for row in _dictfetchall(cursor):
            # Extract fields by name (order-independent, safer than positional unpacking)
            old_id = row["id"]

            measurables = {field: row[field] for field in MEASURABLE_FIELDS}
            measurables["twop_peak_gm"] = measurables.pop("twop_peakGM")

            # Get protein for label
            protein = Protein.objects.get(id=row["protein_id"])

            # Create State (MTI child of Fluorophore)
            state = State.objects.create(
                id=old_id,  # Preserve old ID for easier FK updates later
                created=row["created"],
                modified=row["modified"],
                name=row["name"],
                slug=row["slug"],
                entity_type="p",
                owner_name=protein.name,
                owner_slug=protein.slug,
                created_by_id=row["created_by_id"],
                updated_by_id=row["updated_by_id"],
                **measurables,
                # State-specific fields
                protein_id=protein.id,
                maturation=row["maturation"],
            )

            FluorescenceMeasurement.objects.create(
                id=old_id,  # Preserve old ID
                fluorophore=state,  # State is-a Fluorophore (MTI)
                reference_id=protein.primary_reference_id,
                **measurables,
                is_trusted=True,  # Mark as trusted since it's the original data
                created_by_id=row["created_by_id"],
                updated_by_id=row["updated_by_id"],
            )

    print(f"Migrated {State.objects.count()} State records")

    # Reset the Fluorophore and FluorescenceMeasurement ID sequences
    # to avoid conflicts when creating Dye fluorophores and their measurements
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT setval(
                pg_get_serial_sequence('proteins_fluorophore', 'id'),
                COALESCE((SELECT MAX(id) FROM proteins_fluorophore), 1)
            )
        """)
        print(f"Reset Fluorophore ID sequence to {cursor.fetchone()[0]}")

        cursor.execute("""
            SELECT setval(
                pg_get_serial_sequence('proteins_fluorescencemeasurement', 'id'),
                COALESCE((SELECT MAX(id) FROM proteins_fluorescencemeasurement), 1)
            )
        """)
        print(f"Reset FluorescenceMeasurement ID sequence to {cursor.fetchone()[0]}")


def migrate_dye_data(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> dict[int, int]:
    """Migrate Dye data from old schema to new Dye container + DyeState structure.

    Returns a mapping of old Dye ID → new Fluorophore ID for efficient FK updates.
    """
    Dye = apps.get_model("proteins", "Dye")
    DyeState = apps.get_model("proteins", "DyeState")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    old_to_new_id_map = {}

    # Access old Dye data via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy, emhex, exhex,
                   manufacturer, part, url, created_by_id, updated_by_id
            FROM proteins_dye_old
        """)

        for row in _dictfetchall(cursor):
            old_id = row["id"]
            measurables = {field: row[field] for field in MEASURABLE_FIELDS}
            measurables["twop_peak_gm"] = measurables.pop("twop_peakGM")

            # Create Dye container (without fluorescence properties)
            dye = Dye.objects.create(
                created=row["created"],
                created_by_id=row["created_by_id"],
                updated_by_id=row["updated_by_id"],
                modified=row["modified"],
                name=row["name"],
                slug=row["slug"],
                manufacturer=row["manufacturer"],
                part=row["part"],
                url=row["url"],
            )

            # Create DyeState (MTI child of Fluorophore)
            dyestate = DyeState.objects.create(
                created=row["created"],
                modified=row["modified"],
                created_by_id=row["created_by_id"],
                updated_by_id=row["updated_by_id"],
                name="default",
                slug=f"{dye.slug}-default",
                entity_type="d",
                **measurables,
                owner_name=dye.name,
                owner_slug=dye.slug,
                dye=dye,
            )

            # Create FluorescenceMeasurement from old Dye data
            # For dyes, we don't have a primary_reference concept in old schema
            FluorescenceMeasurement.objects.create(
                fluorophore=dyestate,  # DyeState is-a Fluorophore (MTI)
                reference_id=None,
                **measurables,
                is_trusted=True,
                created_by_id=row["created_by_id"],
                updated_by_id=row["updated_by_id"],
            )

            dye.default_state = dyestate
            dye.save()

            # Store mapping for efficient FK updates
            old_to_new_id_map[old_id] = dyestate.fluorophore_ptr_id

    print(f"Migrated {Dye.objects.count()} Dye records to Dye containers")
    print(f"Created {DyeState.objects.count()} DyeState records")
    return old_to_new_id_map


def update_spectrum_ownership(apps: Apps, schema_editor: BaseDatabaseSchemaEditor, dye_id_map: dict[int, int]) -> None:
    """Update Spectrum foreign keys to point to new Fluorophore records."""
    with schema_editor.connection.cursor() as cursor:
        # Update spectra owned by States - State ID == Fluorophore ID (preserved)
        cursor.execute("""
            UPDATE proteins_spectrum
            SET owner_fluor_id = owner_state_id
            WHERE owner_state_id IS NOT NULL
        """)
        state_count = cursor.rowcount
        print(f"Updated {state_count} spectra owned by States")

        # Update spectra owned by Dyes - need mapping from old Dye ID to new Fluorophore ID
        # Use temp table for the mapping (reuse pattern from OcFluorEff)
        cursor.execute("""
            CREATE TEMP TABLE dye_spectrum_mapping (
                old_id INTEGER PRIMARY KEY,
                new_id INTEGER NOT NULL
            )
        """)

        if dye_id_map:
            values = ", ".join(
                cursor.mogrify("(%s, %s)", (old_id, new_id)).decode() for old_id, new_id in dye_id_map.items()
            )
            cursor.execute(f"INSERT INTO dye_spectrum_mapping (old_id, new_id) VALUES {values}")

        cursor.execute("""
            UPDATE proteins_spectrum s
            SET owner_fluor_id = m.new_id
            FROM dye_spectrum_mapping m
            WHERE s.owner_dye_id = m.old_id
        """)
        dye_count = cursor.rowcount
        print(f"Updated {dye_count} spectra owned by Dyes")

        cursor.execute("DROP TABLE dye_spectrum_mapping")


def update_ocfluoreff(apps: Apps, schema_editor: BaseDatabaseSchemaEditor, dye_id_map: dict[int, int]) -> None:
    """Update OcFluorEff to use direct FK to Fluorophore.

    Uses efficient bulk updates:
    - States: Direct assignment (State ID == Fluorophore ID)
    - Dyes: Temp table JOIN (using provided old_id → new_id mapping)
    """
    with schema_editor.connection.cursor() as cursor:
        # Increase work_mem for efficient hash joins on large tables
        cursor.execute("SET LOCAL work_mem = '256MB'")

        # Cache content_type IDs to avoid repeated subqueries
        cursor.execute("""
            SELECT id FROM django_content_type
            WHERE app_label = 'proteins' AND model = 'state'
        """)
        result = cursor.fetchone()
        if result is None:
            # No content types exist yet (fresh database), nothing to migrate
            return
        state_ct_id = result[0]

        cursor.execute("""
            SELECT id FROM django_content_type
            WHERE app_label = 'proteins' AND model = 'dye'
        """)
        result = cursor.fetchone()
        if result is None:
            return
        dye_ct_id = result[0]

        # First, clean up orphaned records that point to deleted States/Dyes
        # These have object_ids that don't exist in the old tables anymore
        cursor.execute(
            """
            DELETE FROM proteins_ocfluoreff o
            WHERE o.content_type_id = %s
            AND NOT EXISTS (
                SELECT 1 FROM proteins_state_old s WHERE s.id = o.object_id
            )
        """,
            [state_ct_id],
        )
        orphaned_states = cursor.rowcount

        cursor.execute(
            """
            DELETE FROM proteins_ocfluoreff o
            WHERE o.content_type_id = %s
            AND NOT EXISTS (
                SELECT 1 FROM proteins_dye_old d WHERE d.id = o.object_id
            )
        """,
            [dye_ct_id],
        )
        orphaned_dyes = cursor.rowcount
        print(f"Deleted {orphaned_states + orphaned_dyes} orphaned OcFluorEff records")

        # Update OcFluorEff records that pointed to States via GenericFK
        # Since State IDs were preserved, object_id directly maps to fluor_id
        cursor.execute(
            """
            UPDATE proteins_ocfluoreff
            SET fluor_id = object_id
            WHERE content_type_id = %s
        """,
            [state_ct_id],
        )

        state_count = cursor.rowcount
        print(f"Updated {state_count} OcFluorEff records that pointed to States")

        # Update OcFluorEff records that pointed to Dyes via GenericFK
        # Use temp table with primary key for efficient hash join
        cursor.execute("""
            CREATE TEMP TABLE dye_id_mapping (
                old_id INTEGER PRIMARY KEY,
                new_id INTEGER NOT NULL
            )
        """)

        # Bulk insert using execute with VALUES - faster than executemany
        if dye_id_map:
            values = ", ".join(
                cursor.mogrify("(%s, %s)", (old_id, new_id)).decode() for old_id, new_id in dye_id_map.items()
            )
            cursor.execute(f"INSERT INTO dye_id_mapping (old_id, new_id) VALUES {values}")

        cursor.execute(
            """
            UPDATE proteins_ocfluoreff o
            SET fluor_id = m.new_id
            FROM dye_id_mapping m
            WHERE o.object_id = m.old_id
            AND o.content_type_id = %s
        """,
            [dye_ct_id],
        )

        dye_count = cursor.rowcount
        print(f"Updated {dye_count} OcFluorEff records that pointed to Dyes")

        # Clean up temp table
        cursor.execute("DROP TABLE dye_id_mapping")


def populate_emhex_exhex(apps, _schema_editor):
    """Populate emhex/exhex for all Fluorophore objects.

    The historical models don't include the save() logic from AbstractFluorescenceData,
    so we need to manually calculate these fields after creating the objects.
    """

    def wave_to_hex(wavelength, gamma=1):
        """This converts a given wavelength into an approximate RGB value."""
        if not wavelength:
            return "#000"

        wavelength = float(wavelength)
        if 520 <= wavelength:
            wavelength += 40

        if wavelength < 380:
            r = 0.05
            g = 0.0
            b = 0.15
        elif wavelength >= 380 and wavelength <= 440:
            attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
            r = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
            g = 0.0
            b = (1.0 * attenuation) ** gamma
        elif wavelength >= 440 and wavelength <= 490:
            r = 0.0
            g = ((wavelength - 440) / (490 - 440)) ** gamma
            b = 1.0
        elif wavelength >= 490 and wavelength <= 510:
            r = 0.0
            g = 1.0
            b = (-(wavelength - 510) / (510 - 490)) ** gamma
        elif wavelength >= 510 and wavelength <= 580:
            r = ((wavelength - 510) / (580 - 510)) ** gamma
            g = 1.0
            b = 0.0
        elif wavelength >= 580 and wavelength <= 645:
            r = 1.0
            g = (-(wavelength - 645) / (645 - 580)) ** gamma
            b = 0.0
        elif wavelength >= 645 and wavelength <= 750:
            attenuation = 0.3 + 0.7 * (770 - wavelength) / (770 - 645)
            r = (1.0 * attenuation) ** gamma
            g = 0.0
            b = 0.0
        else:
            r = 0.18
            g = 0.0
            b = 0.05
        r *= 255
        g *= 255
        b *= 255
        return f"#{int(r):02x}{int(g):02x}{int(b):02x}"

    Fluorophore = apps.get_model("proteins", "Fluorophore")

    fluorophores_to_update = []
    for fluor in Fluorophore.objects.all():
        fluor.emhex = "#000" if fluor.is_dark else wave_to_hex(fluor.em_max)
        fluor.exhex = wave_to_hex(fluor.ex_max)
        fluorophores_to_update.append(fluor)

    Fluorophore.objects.bulk_update(fluorophores_to_update, ["emhex", "exhex"], batch_size=500)
    print(f"Populated emhex/exhex for {len(fluorophores_to_update)} fluorophores")


def migrate_reversion_state_versions(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Migrate django-reversion State versions to work with new MTI structure.

    Background:
    -----------
    State now inherits from Fluorophore via Multi-Table Inheritance (MTI).
    When reversion reverts a revision, it needs Version records for BOTH the parent
    (Fluorophore) and child (State) models in the same revision.

    Old State versions stored ALL fields in one record:
        {"model": "proteins.state", "fields": {"name": "x", "ex_max": 488, "protein": 1, ...}}

    After migration, we need TWO records per revision:
        1. Fluorophore: {"fields": {"name": "x", "ex_max": 488, ...}}
        2. State: {"fields": {"fluorophore_ptr": 402, "protein": 1, "maturation": 13.0}}

    Both records share the same revision_id, so revision.revert() restores them together.
    """
    ContentType = apps.get_model("contenttypes", "ContentType")

    # Fields that remain on State model; everything else moves to Fluorophore parent
    state_only_fields = {"protein", "maturation", "transitions"}

    with schema_editor.connection.cursor() as cursor:
        try:
            state_ct = ContentType.objects.get(app_label="proteins", model="state")
        except ContentType.DoesNotExist:
            # No content types exist yet (fresh database), nothing to migrate
            return
        fluor_ct, _ = ContentType.objects.get_or_create(app_label="proteins", model="fluorophore")

        # Build lookup dict for protein name/slug (used for owner_name/owner_slug)
        cursor.execute("SELECT id, name, slug FROM proteins_protein")
        proteins = {row[0]: (row[1], row[2]) for row in cursor.fetchall()}

        fluor_inserts: list[tuple[str, str, str, str, str, str, str]] = []
        state_updates: list[tuple[str, str]] = []

        # Fetch all existing State version records
        cursor.execute(
            "SELECT id, object_id, revision_id, serialized_data, db, format "
            "FROM reversion_version WHERE content_type_id = %s",
            [state_ct.id],
        )
        for row in _dictfetchall(cursor):
            version_id = row["id"]
            object_id = row["object_id"]
            revision_id = row["revision_id"]
            serialized_data = row["serialized_data"]
            db = row["db"]
            fmt = row["format"]

            try:
                data = json.loads(serialized_data)
                if not data:
                    continue
                first_record: dict = data[0]
                # data is a list of versioned objects
                old_fields = first_record.get("fields", {})
                pk = first_record.get("pk", object_id)

                # Look up owner (protein) name/slug from the protein FK
                protein_id = old_fields.get("protein")
                owner_name, owner_slug = proteins.get(protein_id, ("", ""))

                # Split fields: State-only fields stay on State, rest go to Fluorophore
                # State needs fluorophore_ptr to link to its MTI parent
                state_fields = {"fluorophore_ptr": pk}
                # Fluorophore needs new required fields with sensible defaults
                fluor_fields = {
                    "entity_type": "p",
                    "owner_name": owner_name,
                    "owner_slug": owner_slug,
                    "source_map": "{}",
                }

                for k, v in old_fields.items():
                    (state_fields if k in state_only_fields else fluor_fields)[k] = v

                # Queue NEW Fluorophore version (same revision_id links them together)
                fluor_inserts.append(
                    (
                        str(object_id),
                        revision_id,  # Same revision as the State version
                        fluor_ct.id,
                        json.dumps([{"model": "proteins.fluorophore", "pk": pk, "fields": fluor_fields}]),
                        "",  # object_repr - not critical for historical versions
                        db,
                        fmt,
                    )
                )
                # Queue UPDATE to existing State version (strip out fields that moved to Fluorophore)
                state_updates.append(
                    (
                        json.dumps([{"model": "proteins.state", "pk": pk, "fields": state_fields}]),
                        version_id,
                    )
                )
            except (json.JSONDecodeError, KeyError, TypeError):
                logger.warning(f"Skipping malformed reversion State version ID {row['id']}")
                continue

        # Bulk insert new Fluorophore versions
        if fluor_inserts:
            cursor.executemany(
                "INSERT INTO reversion_version "
                "(object_id, revision_id, content_type_id, serialized_data, object_repr, db, format) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s)",
                fluor_inserts,
            )
        # Bulk update existing State versions to only contain State-specific fields
        if state_updates:
            cursor.executemany(
                "UPDATE reversion_version SET serialized_data = %s WHERE id = %s",
                state_updates,
            )
        print(f"Migrated {len(fluor_inserts)} State versions for MTI")


def migrate_forward(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Run all migration functions."""
    print("Starting data migration from old schema...")
    migrate_state_data(apps, schema_editor)
    dye_id_map = migrate_dye_data(apps, schema_editor)
    update_spectrum_ownership(apps, schema_editor, dye_id_map)
    update_ocfluoreff(apps, schema_editor, dye_id_map)
    populate_emhex_exhex(apps, schema_editor)
    migrate_reversion_state_versions(apps, schema_editor)
    print("Data migration complete!")


def migrate_reverse(_apps, _schema_editor):
    """Reverse migration - not supported.

    This migration performs a one-way data transformation from the old schema
    (separate State and Dye models) to the new schema (MTI-based Fluorophore hierarchy).

    Reversing this migration would require:
    1. Decomposing Fluorophore + State back into old State structure
    2. Decomposing Fluorophore + DyeState back into old Dye structure
    3. Merging FluorescenceMeasurement data back into parent entities
    4. Restoring old spectrum ownership relationships

    This is not safely automatable. If you need to rollback this migration:
    1. Restore from a database backup taken before running this migration
    2. Do NOT attempt to use Django's reverse migration functionality
    3. Estimated restore time: 30-60 minutes depending on database size

    See deployment documentation for rollback procedures.
    """
    raise RuntimeError(
        "This migration cannot be reversed. "
        "Restore from database backup if rollback is needed. "
        "See migration docstring for details."
    )


def abstract_fluorescence_data_fields():
    """Return fresh field instances for AbstractFluorescenceData.

    Each model needs its own unique field instances to avoid state conflicts.
    """
    return [
        (
            "is_dark",
            models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluoresce"),
        ),
        (
            "ex_max",
            models.PositiveSmallIntegerField(
                blank=True,
                null=True,
                validators=[MinValueValidator(300), MaxValueValidator(900)],
                db_index=True,
                help_text="Excitation maximum (nm)",
            ),
        ),
        (
            "em_max",
            models.PositiveSmallIntegerField(
                blank=True,
                null=True,
                validators=[MinValueValidator(300), MaxValueValidator(1000)],
                db_index=True,
                help_text="Emission maximum (nm)",
            ),
        ),
        (
            "ext_coeff",
            models.IntegerField(
                blank=True,
                null=True,
                verbose_name="Extinction Coefficient (M-1 cm-1)",
                validators=[MinValueValidator(0), MaxValueValidator(300000)],
            ),
        ),
        (
            "qy",
            models.FloatField(
                null=True,
                blank=True,
                verbose_name="Quantum Yield",
                validators=[MinValueValidator(0), MaxValueValidator(1)],
            ),
        ),
        ("brightness", models.FloatField(null=True, blank=True, editable=False)),
        (
            "lifetime",
            models.FloatField(
                null=True,
                blank=True,
                help_text="Lifetime (ns)",
                validators=[MinValueValidator(0), MaxValueValidator(20)],
            ),
        ),
        (
            "pka",
            models.FloatField(
                null=True, blank=True, verbose_name="pKa", validators=[MinValueValidator(2), MaxValueValidator(12)]
            ),
        ),
        (
            "twop_ex_max",
            models.PositiveSmallIntegerField(
                blank=True,
                null=True,
                verbose_name="Peak 2P excitation",
                validators=[MinValueValidator(700), MaxValueValidator(1600)],
                db_index=True,
            ),
        ),
        (
            "twop_peak_gm",
            models.FloatField(
                null=True,
                blank=True,
                verbose_name="Peak 2P cross-section of S0->S1 (GM)",
                validators=[MinValueValidator(0), MaxValueValidator(200)],
            ),
        ),
        (
            "twop_qy",
            models.FloatField(
                null=True,
                blank=True,
                verbose_name="2P Quantum Yield",
                validators=[MinValueValidator(0), MaxValueValidator(1)],
            ),
        ),
        ("emhex", models.CharField(max_length=7, blank=True)),
        ("exhex", models.CharField(max_length=7, blank=True)),
    ]


def authorable_mixin_fields():
    """Return fresh field instances for Authorable mixin."""
    return [
        (
            "created_by",
            models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name="%(class)s_author",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
        (
            "updated_by",
            models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name="%(class)s_modifier",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
    ]


def timestamped_mixin_fields():
    """Return fresh field instances for TimeStampedModel mixin."""
    return [
        (
            "created",
            model_utils.fields.AutoCreatedField(
                default=django.utils.timezone.now, editable=False, verbose_name="created"
            ),
        ),
        (
            "modified",
            model_utils.fields.AutoLastModifiedField(
                default=django.utils.timezone.now, editable=False, verbose_name="modified"
            ),
        ),
    ]


def product_mixin_fields():
    """Return fresh field instances for Product mixin."""
    return [
        ("manufacturer", models.CharField(max_length=128, blank=True)),
        ("part", models.CharField(max_length=128, blank=True)),
        ("url", models.URLField(blank=True)),
    ]


class Migration(migrations.Migration):
    dependencies = [
        ("proteins", "0058_snapgeneplasmid_protein_snapgene_plasmids"),
        ("references", "0001_initial"),
    ]

    operations = [
        # Step 1: move State/Dye tables to _old versions and delete the models.
        migrations.SeparateDatabaseAndState(
            database_operations=[
                # In the database: copy old tables to _old versions, then drop originals
                # This avoids index naming conflicts
                migrations.RunSQL(
                    sql="""
                        CREATE TABLE proteins_state_old AS SELECT * FROM proteins_state;
                        DROP TABLE proteins_state CASCADE;
                    """,
                    reverse_sql="DROP TABLE IF EXISTS proteins_state_old CASCADE;",
                ),
                migrations.RunSQL(
                    sql="""
                        CREATE TABLE proteins_dye_old AS SELECT * FROM proteins_dye;
                        DROP TABLE proteins_dye CASCADE;
                    """,
                    reverse_sql="DROP TABLE IF EXISTS proteins_dye_old CASCADE;",
                ),
            ],
            state_operations=[
                # In Django's state: remove old models
                migrations.DeleteModel(name="State"),
                migrations.DeleteModel(name="Dye"),
            ],
        ),
        # Step 2: Create all new models
        # Fluorophore is the new "State" base model for MTI
        # State and DyeState are MTI children of Fluorophore
        # What was Dye is now is a container model for DyeStates (like Protein for States)
        migrations.CreateModel(
            name="Fluorophore",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                *timestamped_mixin_fields(),
                ("name", models.CharField(max_length=100, db_index=True, default="default")),
                ("slug", models.SlugField(max_length=200, unique=True)),
                (
                    "entity_type",
                    models.CharField(max_length=2, choices=[("p", "Protein"), ("d", "Dye")], db_index=True),
                ),
                (
                    "owner_name",
                    models.CharField(
                        blank=True,
                        db_index=True,
                        default="",
                        help_text="Protein/Dye name (cached for searching)",
                        max_length=255,
                    ),
                ),
                (
                    "owner_slug",
                    models.SlugField(
                        max_length=200, blank=True, default="", help_text="Protein/Dye slug (cached for URLs)"
                    ),
                ),
                *abstract_fluorescence_data_fields(),
                ("source_map", models.JSONField(default=dict, blank=True)),
                ("pinned_source_map", models.JSONField(default=dict, blank=True)),
                *authorable_mixin_fields(),
            ],
            options={
                "indexes": [
                    models.Index(fields=["ex_max"], name="fluorophore_ex_max_idx"),
                    models.Index(fields=["em_max"], name="fluorophore_em_max_idx"),
                    models.Index(fields=["owner_name"], name="fluorophore_owner_name_idx"),
                    models.Index(fields=["entity_type", "is_dark"], name="fluorophore_type_dark_idx"),
                ],
            },
        ),
        migrations.CreateModel(
            name="FluorescenceMeasurement",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                *abstract_fluorescence_data_fields(),
                ("date_measured", models.DateField(null=True, blank=True)),
                ("conditions", models.TextField(blank=True, help_text="pH, solvent, temp, etc.")),
                (
                    "is_trusted",
                    models.BooleanField(default=False, help_text="If True, this measurement overrides others."),
                ),
                (
                    "fluorophore",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="measurements",
                        to="proteins.fluorophore",
                    ),
                ),
                (
                    "reference",
                    models.ForeignKey(
                        "references.Reference", on_delete=django.db.models.deletion.CASCADE, null=True, blank=True
                    ),
                ),
                *authorable_mixin_fields(),
                *timestamped_mixin_fields(),
            ],
            options={
                "abstract": False,
            },
        ),
        migrations.CreateModel(
            name="Dye",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("name", models.CharField(max_length=255, db_index=True)),
                ("slug", models.SlugField(unique=True)),
                (
                    "primary_reference",
                    models.ForeignKey(
                        blank=True,
                        help_text="The publication that introduced the dye",
                        null=True,
                        on_delete=django.db.models.deletion.SET_NULL,
                        related_name="primary_dyes",
                        to="references.reference",
                        verbose_name="Primary Reference",
                    ),
                ),
                *product_mixin_fields(),
                *authorable_mixin_fields(),
                *timestamped_mixin_fields(),
            ],
        ),
        migrations.CreateModel(
            name="DyeState",
            fields=[
                (
                    "fluorophore_ptr",
                    models.OneToOneField(
                        auto_created=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        parent_link=True,
                        primary_key=True,
                        serialize=False,
                        to="proteins.fluorophore",
                    ),
                ),
                (
                    "dye",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE, related_name="states", to="proteins.dye"
                    ),
                ),
            ],
            options={
                "abstract": False,
            },
            bases=("proteins.fluorophore",),
        ),
        migrations.AddField(
            model_name="dye",
            name="default_state",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name="default_for",
                to="proteins.DyeState",
            ),
        ),
        # Re-create State model as MTI child of Fluorophore
        migrations.CreateModel(
            name="State",
            fields=[
                (
                    "fluorophore_ptr",
                    models.OneToOneField(
                        auto_created=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        parent_link=True,
                        primary_key=True,
                        serialize=False,
                        to="proteins.fluorophore",
                    ),
                ),
                (
                    "protein",
                    models.ForeignKey(
                        help_text="The protein to which this state belongs",
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="states",
                        to="proteins.protein",
                    ),
                ),
                (
                    "maturation",
                    models.FloatField(
                        null=True,
                        blank=True,
                        help_text="Maturation time (min)",
                        validators=[MinValueValidator(0), MaxValueValidator(1600)],
                    ),
                ),
                (
                    "transitions",
                    models.ManyToManyField(
                        blank=True,
                        related_name="transition_state",
                        through="proteins.StateTransition",
                        to="proteins.state",
                        verbose_name="State Transitions",
                    ),
                ),
            ],
            options={
                "abstract": False,
            },
            bases=("proteins.fluorophore",),
        ),
        # Add owner_fluor to Spectrum (nullable for now)
        migrations.AddField(
            model_name="spectrum",
            name="owner_fluor",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="spectra",
                to="proteins.fluorophore",
            ),
        ),
        # Add fluor FK to OcFluorEff (nullable for now)
        migrations.AddField(
            model_name="ocfluoreff",
            name="fluor",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="oc_effs",
                to="proteins.fluorophore",
            ),
        ),
        # Add index for spectrum owner_fluor lookups
        migrations.AddIndex(
            model_name="spectrum",
            index=models.Index(fields=["owner_fluor_id", "status"], name="spectrum_fluor_status_idx"),
        ),
        # ----------------------------------------------
        # Perform manual data migration steps
        migrations.RunPython(migrate_forward, migrate_reverse),
        # ----------------------------------------------
        # Step 1: Drop old tables that we renamed in 0059
        migrations.RunSQL(
            sql="DROP TABLE IF EXISTS proteins_state_old CASCADE;",
            reverse_sql=migrations.RunSQL.noop,
        ),
        migrations.RunSQL(
            sql="DROP TABLE IF EXISTS proteins_dye_old CASCADE;",
            reverse_sql=migrations.RunSQL.noop,
        ),
        # Step 2: Remove old foreign keys from Spectrum
        migrations.RemoveField(model_name="spectrum", name="owner_state"),
        migrations.RemoveField(model_name="spectrum", name="owner_dye"),
        migrations.RemoveIndex(model_name="spectrum", name="spectrum_state_status_idx"),
        # Step 3: Make owner_fluor non-nullable now that all data is migrated
        # First, verify no nulls exist (will fail if there are any)
        migrations.RunSQL(
            sql="""
                DO $$
                BEGIN
                    IF EXISTS (
                        SELECT 1 FROM proteins_spectrum
                        WHERE owner_fluor_id IS NULL
                        AND (owner_filter_id IS NULL
                             AND owner_light_id IS NULL
                             AND owner_camera_id IS NULL)
                    ) THEN
                        RAISE EXCEPTION 'Found Spectrum records with no owner after migration!';
                    END IF;
                END $$;
            """,
            reverse_sql=migrations.RunSQL.noop,
        ),
        # Step 4: Remove GenericForeignKey fields from OcFluorEff
        migrations.RemoveField(
            model_name="ocfluoreff",
            name="content_type",
        ),
        migrations.RemoveField(
            model_name="ocfluoreff",
            name="object_id",
        ),
        # Step 5: Make fluor FK on OcFluorEff non-nullable
        migrations.AlterField(
            model_name="ocfluoreff",
            name="fluor",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="oc_effs",
                to="proteins.fluorophore",
            ),
        ),
        migrations.AlterUniqueTogether(name="ocfluoreff", unique_together={("oc", "fluor")}),
    ]
