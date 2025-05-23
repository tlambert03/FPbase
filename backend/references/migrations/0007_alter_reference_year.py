# Generated by Django 4.2.1 on 2024-01-04 16:39

import django.core.validators
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("references", "0006_alter_author_id_alter_reference_id_and_more"),
    ]

    operations = [
        migrations.AlterField(
            model_name="reference",
            name="year",
            field=models.PositiveIntegerField(
                help_text="YYYY",
                validators=[
                    django.core.validators.MinLengthValidator(4),
                    django.core.validators.MaxLengthValidator(4),
                    django.core.validators.MinValueValidator(1960),
                    django.core.validators.MaxValueValidator(2025),
                ],
            ),
        ),
    ]
