# Generated by Django 2.0.5 on 2018-05-12 00:58

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('references', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='reference',
            name='firstauthor',
            field=models.CharField(blank=True, default='', max_length=100),
        ),
    ]