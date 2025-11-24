from __future__ import annotations
from typing import Any
# Generated manually for schema overhaul

from django.core.validators import MaxValueValidator, MinValueValidator
from django.contrib.postgres.fields import ArrayField
from django.db import migrations, models
import django.db.models.deletion
import model_utils.fields
import django.utils.timezone
import logging
from django.db import migrations, models
import django.db.models.deletion
from django.db import migrations
from django.apps.registry import Apps
from django.db.backends.base.schema import BaseDatabaseSchemaEditor
from django.db.backends.utils import CursorWrapper
logger = logging.getLogger(__name__)


def _dictfetchall(cursor: CursorWrapper) -> list[dict[str, Any]]:
    """Return all rows from a cursor as a dict. Assume the column names are unique."""
    if not cursor.description:
        return []
    columns = (col.name for col in cursor.description)
    return [dict(zip(columns, row)) for row in cursor.fetchall()]

def migrate_state_data(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Migrate State data from old schema to new Fluorophore + State MTI structure."""
    # Get models from migration state
    Fluorophore = apps.get_model("proteins", "Fluorophore")
    State = apps.get_model("proteins", "State")
    Protein = apps.get_model("proteins", "Protein")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    # Access old State data directly via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy,
                   maturation, protein_id, created_by_id, updated_by_id
            FROM proteins_state_old
        """)

        for row in cursor.fetchall():
            (old_id, created, modified, name, slug, is_dark,
             ex_max, em_max, ext_coeff, qy, brightness,
             lifetime, pka, twop_ex_max, twop_peak_gm, twop_qy,
             maturation, protein_id, created_by_id, updated_by_id) = row

            # Get protein for label
            try:
                protein = Protein.objects.get(id=protein_id)
            except Protein.DoesNotExist:
                # note, the db doesn't have any cases of this
                print(f"Warning: Protein {protein_id} not found for State {old_id}, skipping")
                continue

            # Handle empty/null slugs with guaranteed non-empty fallback
            if not slug or (isinstance(slug, str) and slug.strip() == ""):
                # Try protein slug + name, or protein slug, or fallback to state ID
                if name and name != "default":
                    state_slug = f"{protein.slug}-{name}" if protein.slug else f"state-{old_id}"
                else:
                    state_slug = protein.slug if protein.slug else f"state-{old_id}"
            else:
                state_slug = slug

            # Final safety check - ensure slug is not empty
            if not state_slug or state_slug.strip() == "":
                state_slug = f"state-{old_id}"


            # Ensure slug uniqueness by checking if it already exists
            original_slug = state_slug
            base_slug = state_slug
            counter = 1
            while Fluorophore.objects.filter(slug=state_slug).exists():
                state_slug = f"{base_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique slug for State {old_id} "
                        f"after 100 attempts (original: {original_slug})"
                    )

            if state_slug != original_slug:
                logger.warning(
                    f"Slug collision during State migration: {original_slug} -> {state_slug} "
                    f"(State ID: {old_id}, Protein: {protein.name})"
                )

            # Create Fluorophore parent (MTI will link automatically)
            fluorophore = Fluorophore.objects.create(
                id=old_id,  # Preserve old ID for easier FK updates later
                created=created,
                modified=modified,
                name=name,
                slug=state_slug,
                entity_type="p",
                owner_name=protein.name,
                owner_slug=protein.slug,
                ex_max=ex_max,
                em_max=em_max,
                ext_coeff=ext_coeff,
                qy=qy,
                brightness=brightness,
                lifetime=lifetime,
                pka=pka,
                twop_ex_max=twop_ex_max,
                twop_peak_gm=twop_peak_gm,  # Map from SQL result to model field
                twop_qy=twop_qy,
                is_dark=is_dark,
                created_by_id=created_by_id,
                updated_by_id=updated_by_id
            )

            # Create State (MTI child) pointing to the Fluorophore
            # Use raw SQL to avoid Django MTI trying to update parent with empty values
            cursor.execute("""
                INSERT INTO proteins_state (fluorophore_ptr_id, protein_id, maturation)
                VALUES (%s, %s, %s)
            """, [fluorophore.pk, protein_id, maturation])

            # Create FluorescenceMeasurement from old State data if there's any fluorescence data
            if any([ex_max, em_max, qy, ext_coeff, lifetime, pka, twop_ex_max, twop_peak_gm, twop_qy]):
                FluorescenceMeasurement.objects.create(
                    id=old_id,  # Preserve old ID
                    fluorophore=fluorophore,
                    reference_id=protein.primary_reference_id,
                    ex_max=ex_max,
                    em_max=em_max,
                    ext_coeff=ext_coeff,
                    qy=qy,
                    brightness=brightness,
                    lifetime=lifetime,
                    pka=pka,
                    twop_ex_max=twop_ex_max,
                    twop_peak_gm=twop_peak_gm,
                    twop_qy=twop_qy,
                    is_dark=is_dark,
                    is_trusted=True,  # Mark as trusted since it's the original data
                    created_by_id=created_by_id,
                    updated_by_id=updated_by_id
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
        fluor_seq = cursor.fetchone()[0]
        print(f"Reset Fluorophore ID sequence to {fluor_seq}")

        cursor.execute("""
            SELECT setval(
                pg_get_serial_sequence('proteins_fluorescencemeasurement', 'id'),
                COALESCE((SELECT MAX(id) FROM proteins_fluorescencemeasurement), 1)
            )
        """)
        meas_seq = cursor.fetchone()[0]
        print(f"Reset FluorescenceMeasurement ID sequence to {meas_seq}")


def migrate_dye_data(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Migrate Dye data from old schema to new Dye container + DyeState structure."""
    Fluorophore = apps.get_model("proteins", "Fluorophore")
    Dye = apps.get_model("proteins", "Dye")
    DyeState = apps.get_model("proteins", "DyeState")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    # Access old Dye data via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy
            FROM proteins_dye_old
        """)

        for row in cursor.fetchall():
            (old_id, created, modified, name, slug, is_dark,
             ex_max, em_max, ext_coeff, qy, brightness,
             lifetime, pka, twop_ex_max, twop_peak_gm, twop_qy) = row

            # Handle empty/null slugs
            if not slug or slug.strip() == "":
                dye_slug = f"dye-{old_id}"  # Use old ID as fallback
            else:
                dye_slug = slug

            # Ensure Dye slug uniqueness
            original_dye_slug = dye_slug
            base_dye_slug = dye_slug
            counter = 1
            while Dye.objects.filter(slug=dye_slug).exists():
                dye_slug = f"{base_dye_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique slug for Dye {old_id} "
                        f"after 100 attempts (original: {original_dye_slug})"
                    )

            if dye_slug != original_dye_slug:
                logger.warning(
                    f"Slug collision during Dye migration: {original_dye_slug} -> {dye_slug} "
                    f"(Dye ID: {old_id}, Name: {name})"
                )

            # Old Dye schema doesn't have chemical structure fields
            # Mark all as PROPRIETARY to avoid unique constraint issues
            # (Can be updated later with actual chemical data)

            # Create Dye container (without fluorescence properties)
            dye = Dye.objects.create(
                created=created,
                modified=modified,
                name=name,
                slug=dye_slug,
            )

            # Create Fluorophore for this DyeState
            # Ensure unique fluorophore slug
            fluorophore_slug = f"{dye_slug}-default"
            original_fluor_slug = fluorophore_slug
            base_fluor_slug = fluorophore_slug
            counter = 1
            while Fluorophore.objects.filter(slug=fluorophore_slug).exists():
                fluorophore_slug = f"{base_fluor_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique fluorophore slug for Dye {old_id} "
                        f"after 100 attempts (original: {original_fluor_slug})"
                    )

            if fluorophore_slug != original_fluor_slug:
                logger.warning(
                    f"Fluorophore slug collision during Dye migration: {original_fluor_slug} -> {fluorophore_slug} "
                    f"(Dye ID: {old_id}, Name: {name})"
                )

            # Don't pass emhex/exhex - the save() method will compute them from wavelengths
            fluorophore = Fluorophore.objects.create(
                created=created,
                modified=modified,
                name="default",
                slug=fluorophore_slug,
                entity_type="d",
                owner_name=dye.name,
                owner_slug=dye.slug,
                ex_max=ex_max,
                em_max=em_max,
                ext_coeff=ext_coeff,
                qy=qy,
                brightness=brightness,
                lifetime=lifetime,
                pka=pka,
                twop_ex_max=twop_ex_max,
                twop_peak_gm=twop_peak_gm,
                twop_qy=twop_qy,
                is_dark=is_dark,
            )

            # Create DyeState (one per old Dye)
            # Use raw SQL to avoid Django MTI trying to update parent with empty values
            cursor.execute("""
                INSERT INTO proteins_dyestate (fluorophore_ptr_id, dye_id)
                VALUES (%s, %s)
            """, [fluorophore.pk, dye.pk])

            # Create FluorescenceMeasurement from old Dye data
            if any([ex_max, em_max, qy, ext_coeff, lifetime, pka]):
                # For dyes, we don't have a primary_reference concept in old schema
                # We'll leave reference as None for now
                FluorescenceMeasurement.objects.create(
                    fluorophore=fluorophore,
                    reference_id=None,
                    ex_max=ex_max,
                    em_max=em_max,
                    ext_coeff=ext_coeff,
                    qy=qy,
                    brightness=brightness,
                    lifetime=lifetime,
                    pka=pka,
                    twop_ex_max=twop_ex_max,
                    twop_peak_gm=twop_peak_gm,
                    twop_qy=twop_qy,
                    is_dark=is_dark,
                    is_trusted=True,
                )

    print(f"Migrated {Dye.objects.count()} Dye records to Dye containers")
    print(f"Created {DyeState.objects.count()} DyeState records")


def update_spectrum_ownership(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Update Spectrum foreign keys to point to new Fluorophore records."""
    Spectrum = apps.get_model("proteins", "Spectrum")
    Fluorophore = apps.get_model("proteins", "Fluorophore")

    # We need to map old State/Dye IDs to new Fluorophore IDs
    # This is tricky because we need to query the old tables

    with schema_editor.connection.cursor() as cursor:
        # Update spectra that were owned by States
        cursor.execute("""
            UPDATE proteins_spectrum s
            SET owner_fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_state ns ON ns.fluorophore_ptr_id = f.id
                JOIN proteins_state_old os ON os.slug = f.slug
                WHERE os.id = s.owner_state_id
            )
            WHERE s.owner_state_id IS NOT NULL
        """)

        state_count = cursor.rowcount
        print(f"Updated {state_count} spectra owned by States")

        # Update spectra that were owned by Dyes
        # Note: DyeState slug is "{dye_slug}-default", so we need to match carefully
        cursor.execute("""
            UPDATE proteins_spectrum s
            SET owner_fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_dyestate ds ON ds.fluorophore_ptr_id = f.id
                JOIN proteins_dye d ON d.id = ds.dye_id
                JOIN proteins_dye_old od ON od.slug = d.slug
                WHERE od.id = s.owner_dye_id
            )
            WHERE s.owner_dye_id IS NOT NULL
        """)

        dye_count = cursor.rowcount
        print(f"Updated {dye_count} spectra owned by Dyes")


def update_ocfluoreff(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Update OcFluorEff to use direct FK to Fluorophore."""
    with schema_editor.connection.cursor() as cursor:
        # Update OcFluorEff records that pointed to States via GenericFK
        cursor.execute("""
            UPDATE proteins_ocfluoreff o
            SET fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_state s ON s.fluorophore_ptr_id = f.id
                WHERE s.fluorophore_ptr_id IN (
                    SELECT fluorophore_ptr_id
                    FROM proteins_state ps
                    JOIN proteins_state_old os ON os.slug = (
                        SELECT slug FROM proteins_fluorophore WHERE id = ps.fluorophore_ptr_id
                    )
                    WHERE os.id = o.object_id::integer
                )
            )
            WHERE o.content_type_id = (
                SELECT id FROM django_content_type
                WHERE app_label = 'proteins' AND model = 'state'
            )
        """)

        state_count = cursor.rowcount
        print(f"Updated {state_count} OcFluorEff records that pointed to States")

        # Update OcFluorEff records that pointed to Dyes via GenericFK
        cursor.execute("""
            UPDATE proteins_ocfluoreff o
            SET fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_dyestate ds ON ds.fluorophore_ptr_id = f.id
                JOIN proteins_dye d ON d.id = ds.dye_id
                JOIN proteins_dye_old od ON od.slug = d.slug
                WHERE od.id = o.object_id::integer
            )
            WHERE o.content_type_id = (
                SELECT id FROM django_content_type
                WHERE app_label = 'proteins' AND model = 'dye'
            )
        """)

        dye_count = cursor.rowcount
        print(f"Updated {dye_count} OcFluorEff records that pointed to Dyes")


def populate_emhex_exhex(apps, _schema_editor):
    """Populate emhex/exhex for all Fluorophore objects.

    The historical models don't include the save() logic from AbstractFluorescenceData,
    so we need to manually calculate these fields after creating the objects.
    """
    from proteins.util.helpers import wave_to_hex

    Fluorophore = apps.get_model("proteins", "Fluorophore")

    fluorophores_to_update = []
    for fluor in Fluorophore.objects.all():
        fluor.emhex = "#000" if fluor.is_dark else wave_to_hex(fluor.em_max)
        fluor.exhex = wave_to_hex(fluor.ex_max)
        fluorophores_to_update.append(fluor)

    Fluorophore.objects.bulk_update(
        fluorophores_to_update,
        ["emhex", "exhex"],
        batch_size=500
    )
    print(f"Populated emhex/exhex for {len(fluorophores_to_update)} fluorophores")


def migrate_forward(apps: Apps, schema_editor: BaseDatabaseSchemaEditor) -> None:
    """Run all migration functions."""
    print("Starting data migration from old schema...")
    migrate_state_data(apps, schema_editor)
    migrate_dye_data(apps, schema_editor)
    update_spectrum_ownership(apps, schema_editor)
    update_ocfluoreff(apps, schema_editor)
    populate_emhex_exhex(apps, schema_editor)
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
        ("is_dark", models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluorescence")),
        ("ex_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(900)], db_index=True, help_text="Excitation maximum (nm)")),
        ("em_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(1000)], db_index=True, help_text="Emission maximum (nm)")),
        ("ext_coeff", models.IntegerField(blank=True, null=True, verbose_name="Extinction Coefficient (M-1 cm-1)", validators=[MinValueValidator(0), MaxValueValidator(300000)])),
        ("qy", models.FloatField(null=True, blank=True, verbose_name="Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
        ("brightness", models.FloatField(null=True, blank=True, editable=False)),
        ("lifetime", models.FloatField(null=True, blank=True, help_text="Lifetime (ns)", validators=[MinValueValidator(0), MaxValueValidator(20)])),
        ("pka", models.FloatField(null=True, blank=True, verbose_name="pKa", validators=[MinValueValidator(2), MaxValueValidator(12)])),
        ("twop_ex_max", models.PositiveSmallIntegerField(blank=True, null=True, verbose_name="Peak 2P excitation", validators=[MinValueValidator(700), MaxValueValidator(1600)], db_index=True)),
        ("twop_peak_gm", models.FloatField(null=True, blank=True, verbose_name="Peak 2P cross-section of S0->S1 (GM)", validators=[MinValueValidator(0), MaxValueValidator(200)])),
        ("twop_qy", models.FloatField(null=True, blank=True, verbose_name="2P Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
        ("emhex", models.CharField(max_length=7, blank=True)),
        ("exhex", models.CharField(max_length=7, blank=True)),
    ]


def authorable_mixin_fields():
    """Return fresh field instances for Authorable mixin."""
    return [
        ("created_by_id", models.IntegerField(null=True, blank=True)),
        ("updated_by_id", models.IntegerField(null=True, blank=True)),
    ]


def timestamped_mixin_fields():
    """Return fresh field instances for TimeStampedModel mixin."""
    return [
        ("created", model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name="created")),
        ("modified", model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name="modified")),
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
                ("entity_type", models.CharField(max_length=2, choices=[("p", "Protein"), ("d", "Dye")], db_index=True)),
                ("owner_name", models.CharField(max_length=255, db_index=True, blank=True, null=True, help_text="Protein/Dye name (cached for searching)")),
                ("owner_slug", models.SlugField(max_length=200, blank=True, null=True, help_text="Protein/Dye slug (cached for URLs)")),
                *abstract_fluorescence_data_fields(),
                ("source_map", models.JSONField(default=dict, blank=True)),
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
                ("is_trusted", models.BooleanField(default=False, help_text="If True, this measurement overrides others.")),
                ("fluorophore", models.ForeignKey("Fluorophore", related_name="measurements", on_delete=django.db.models.deletion.CASCADE)),
                ("reference", models.ForeignKey("references.Reference", on_delete=django.db.models.deletion.CASCADE, null=True, blank=True)),
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
                ("slug", models.SlugField(max_length=100, unique=True)),  # Increased from default 50 to accommodate long names
                *product_mixin_fields(),
                *authorable_mixin_fields(),
                *timestamped_mixin_fields(),
            ],
        ),

        migrations.CreateModel(
            name="DyeState",
            fields=[
                ("fluorophore_ptr", models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to="proteins.fluorophore")),
                ("dye", models.ForeignKey("Dye", on_delete=django.db.models.deletion.CASCADE, related_name="states")),
            ],
            options={
                "abstract": False,
            },
            bases=("proteins.fluorophore",),
        ),

        # Re-create State model as MTI child of Fluorophore
        migrations.CreateModel(
            name="State",
            fields=[
                ("fluorophore_ptr", models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to="proteins.fluorophore")),
                ("protein", models.ForeignKey("Protein", related_name="states", help_text="The protein to which this state belongs", on_delete=django.db.models.deletion.CASCADE)),
                ("maturation", models.FloatField(null=True, blank=True, help_text="Maturation time (min)", validators=[MinValueValidator(0), MaxValueValidator(1600)])),
                ("transitions", models.ManyToManyField(blank=True, related_name="transition_state", through="proteins.StateTransition", to="proteins.state", verbose_name="State Transitions")),
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
        migrations.RemoveField(
            model_name="spectrum",
            name="owner_state",
        ),
        migrations.RemoveField(
            model_name="spectrum",
            name="owner_dye",
        ),

        # Step 3: Make owner_fluor non-nullable now that all data is migrated
        # First, verify no nulls exist (will fail if there are any)
        migrations.RunSQL(
            sql="""
                DO $$
                BEGIN
                    IF EXISTS (SELECT 1 FROM proteins_spectrum WHERE owner_fluor_id IS NULL AND (owner_filter_id IS NULL AND owner_light_id IS NULL AND owner_camera_id IS NULL)) THEN
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
    ]
