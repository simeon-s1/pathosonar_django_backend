# Generated by Django 4.2.4 on 2023-09-19 09:00

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("rest_api", "0006_alignment2mutation_mutation_mutation2annotation_and_more"),
    ]

    operations = [
        migrations.CreateModel(
            name="EnteredData",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("type", models.CharField(blank=True, max_length=50, null=True)),
                ("name", models.CharField(blank=True, max_length=50, null=True)),
                ("date", models.DateField(blank=True, null=True)),
            ],
        ),
    ]