# Generated by Django 2.1.2 on 2018-11-07 20:15

import django.contrib.postgres.fields
from django.db import migrations, models
import proteins.models.lineage
import proteins.validators


class Migration(migrations.Migration):

    dependencies = [
        ('proteins', '0031_auto_20181103_1531'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lineage',
            name='mutation',
            field=proteins.models.lineage.MutationSetField(blank=True, max_length=400),
        ),
        migrations.AlterField(
            model_name='mutation',
            name='mutations',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=5), size=None, validators=[proteins.validators.validate_mutation]),
        ),
    ]
