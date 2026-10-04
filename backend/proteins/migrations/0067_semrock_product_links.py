from django.db import migrations
from django.db.models import F, Value
from django.db.models.functions import Replace

# semrock.com now redirects every FilterDetails page to a 404 on idex-hs.com
OLD = "https://www.semrock.com/FilterDetails.aspx?id="
NEW = "https://www.idex-hs.com/store/search-results/1/?searchCriteria="


def _rewrite(apps, old: str, new: str) -> None:
    # update(), not save(): saving a filter recomputes its slug
    Filter = apps.get_model("proteins", "Filter")
    Filter.objects.filter(url__startswith=old).update(
        url=Replace(F("url"), Value(old), Value(new))
    )


def forward(apps, schema_editor):
    _rewrite(apps, OLD, NEW)


def backward(apps, schema_editor):
    _rewrite(apps, NEW, OLD)


class Migration(migrations.Migration):
    dependencies = [("proteins", "0066_organism_photos")]

    operations = [migrations.RunPython(forward, backward)]
