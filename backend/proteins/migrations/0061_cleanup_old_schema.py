# Generated manually for schema overhaul
#
# WARNING: This migration is NOT REVERSIBLE
# ==========================================
# This migration drops the old State and Dye tables and removes deprecated fields.
# Once this migration runs, there is no automated way to reverse it.
#
# Rollback procedure:
# 1. Restore from database backup taken before migration 0059
# 2. Do NOT use Django's migrate command to reverse - it will not work
# 3. Manual data recovery may be required if backup is unavailable
#
# Always take a full database backup before running this migration in production.

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [
        ("proteins", "0060_migrate_data_from_old_schema"),
    ]

    operations = [
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
