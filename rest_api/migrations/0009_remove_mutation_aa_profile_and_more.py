# Generated by Django 4.2.4 on 2023-09-19 09:45

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("rest_api", "0008_alter_entereddata_name"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="mutation",
            name="aa_profile",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="collection_date",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="country",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="genome_completeness",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="geo_location",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="host",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="imported",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="isolate",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="length",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="nuc_profile",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="reference_accession",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="release_date",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="sample_name",
        ),
        migrations.RemoveField(
            model_name="mutation",
            name="seq_tech",
        ),
    ]
