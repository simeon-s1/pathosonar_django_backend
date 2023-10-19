# Generated by Django 4.2.4 on 2023-09-27 12:11

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("rest_api", "0013_alter_entereddata_unique_together_and_more"),
    ]

    operations = [
        migrations.AddField(
            model_name="sample",
            name="collection_date",
            field=models.DateField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="country",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="genome_completeness",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="host",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="lab",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="length",
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="lineage",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="processing_date",
            field=models.DateField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="sequencing_tech",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="technology",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sample",
            name="zip_code",
            field=models.CharField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name="molecule",
            name="standard",
            field=models.BooleanField(blank=True, null=True),
        ),
    ]
