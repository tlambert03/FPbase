from django.db import migrations


def delete_spectrumless_owners(apps, schema_editor):
    # filters, lights and cameras whose spectrum was deleted (FPBASE-5E6); deleting a
    # spectrum now deletes its owner, too
    for model_name in ("Filter", "Light", "Camera"):
        apps.get_model("proteins", model_name).objects.filter(spectrum__isnull=True).delete()


class Migration(migrations.Migration):
    dependencies = [("proteins", "0061_spectrum_unique_approved_fluor_subtype")]

    operations = [migrations.RunPython(delete_spectrumless_owners, migrations.RunPython.noop)]
