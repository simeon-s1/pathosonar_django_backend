# Generated by Django 4.2.4 on 2023-10-05 10:27

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("rest_api", "0016_alter_gene_table"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="gene",
            name="symbol",
        ),
    ]