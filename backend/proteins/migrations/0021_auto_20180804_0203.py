# Generated by Django 2.0.7 on 2018-08-04 02:03

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('proteins', '0020_auto_20180729_0234'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='fretpair',
            unique_together=set(),
        ),
        migrations.RemoveField(
            model_name='fretpair',
            name='acceptor',
        ),
        migrations.RemoveField(
            model_name='fretpair',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='fretpair',
            name='donor',
        ),
        migrations.RemoveField(
            model_name='fretpair',
            name='pair_references',
        ),
        migrations.RemoveField(
            model_name='fretpair',
            name='updated_by',
        ),
        migrations.RemoveField(
            model_name='protein',
            name='FRET_partner',
        ),
        migrations.DeleteModel(
            name='FRETpair',
        ),
    ]
