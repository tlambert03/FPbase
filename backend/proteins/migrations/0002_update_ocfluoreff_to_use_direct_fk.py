from django.db import migrations, models
import django.db.models.deletion


def migrate_generic_fk_to_direct_fk(apps, schema_editor):
    """Migrate data from GenericForeignKey (content_type/object_id) to direct ForeignKey (fluor)."""
    OcFluorEff = apps.get_model("proteins", "OcFluorEff")
    ContentType = apps.get_model("contenttypes", "ContentType")

    # Get content types for State and DyeState
    try:
        state_ct = ContentType.objects.get(app_label="proteins", model="state")
        dyestate_ct = ContentType.objects.get(app_label="proteins", model="dyestate")
    except ContentType.DoesNotExist:
        # If content types don't exist, there's no data to migrate
        return

    # Migrate all OcFluorEff records
    for eff in OcFluorEff.objects.all():
        # The object_id already points to the correct Fluorophore ID
        # because State and DyeState inherit from Fluorophore via MTI
        if eff.content_type_id in (state_ct.id, dyestate_ct.id):
            eff.fluor_id = eff.object_id
            eff.save(update_fields=['fluor_id'])


class Migration(migrations.Migration):
    dependencies = [
        ("proteins", "0001_initial"),
        ("contenttypes", "0002_remove_content_type_name"),
    ]

    operations = [
        # Step 1: Remove old unique_together constraint first (before removing fields)
        migrations.AlterUniqueTogether(
            name="ocfluoreff",
            unique_together=set(),
        ),
        # Step 2: Add fluor field as nullable
        migrations.AddField(
            model_name="ocfluoreff",
            name="fluor",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                to="proteins.fluorophore",
            ),
        ),
        # Step 3: Migrate data
        migrations.RunPython(
            migrate_generic_fk_to_direct_fk,
            reverse_code=migrations.RunPython.noop,
        ),
        # Step 4: Make fluor non-nullable
        migrations.AlterField(
            model_name="ocfluoreff",
            name="fluor",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                to="proteins.fluorophore",
            ),
        ),
        # Step 5: Remove content_type and object_id fields
        migrations.RemoveField(
            model_name="ocfluoreff",
            name="content_type",
        ),
        migrations.RemoveField(
            model_name="ocfluoreff",
            name="object_id",
        ),
        # Step 6: Add new unique_together constraint
        migrations.AlterUniqueTogether(
            name="ocfluoreff",
            unique_together={("oc", "fluor")},
        ),
    ]
