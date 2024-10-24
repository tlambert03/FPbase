# Generated by Django 2.1.2 on 2019-01-21 13:41

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('proteins', '0045_dye_is_dark'),
    ]

    operations = [
        migrations.AlterField(
            model_name='protein',
            name='pdb',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=4), blank=True, null=True, size=None, verbose_name='Protein DataBank IDs'),
        ),
        migrations.AlterField(
            model_name='protein',
            name='switch_type',
            field=models.CharField(blank=True, choices=[('b', 'Basic'), ('pa', 'Photoactivatable'), ('ps', 'Photoswitchable'), ('pc', 'Photoconvertible'), ('mp', 'Multi-photochromic'), ('o', 'Multistate'), ('t', 'Timer')], default='b', help_text='Photoswitching type (basic if none)', max_length=2, verbose_name='Switching Type'),
        ),
    ]
