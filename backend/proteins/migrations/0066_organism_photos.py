import json
from pathlib import Path

from django.db import migrations

# Curated Wikimedia Commons photos; each was checked against its Commons description
DATA = Path(__file__).with_suffix(".json")


def add_photos(apps, schema_editor):
    Organism = apps.get_model("proteins", "Organism")
    OrganismPhoto = apps.get_model("proteins", "OrganismPhoto")
    photos = json.loads(DATA.read_text())
    existing = set(Organism.objects.values_list("id", flat=True))
    OrganismPhoto.objects.bulk_create(
        [OrganismPhoto(**p) for p in photos if p["organism_id"] in existing],
        ignore_conflicts=True,
    )


def remove_photos(apps, schema_editor):
    apps.get_model("proteins", "OrganismPhoto").objects.all().delete()


class Migration(migrations.Migration):
    dependencies = [("proteins", "0065_organismphoto")]

    operations = [migrations.RunPython(add_photos, remove_photos)]
