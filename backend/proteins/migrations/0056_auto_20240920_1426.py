# Generated by Django 4.2.1 on 2024-09-20 14:26
import logging
from contextlib import suppress
from django.db import migrations
from proteins.models import Spectrum

logger = logging.getLogger(__name__)


def set_default_status(apps, schema_editor):
    for spectrum in Spectrum.objects.all_objects():
        try:
            spectrum.status = Spectrum.STATUS.approved
            spectrum.save()
        except Exception as e:
            logger.error(f"Failed to resave spectrum id: {spectrum.id}", exc_info=e)

class Migration(migrations.Migration):

    dependencies = [
        ("proteins", "0055_spectrum_status_spectrum_status_changed"),
    ]

    operations = [migrations.RunPython(set_default_status)]