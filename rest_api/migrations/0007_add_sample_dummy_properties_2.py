import random
from django.db import migrations
from rest_api.models import Sample, Property, Sample2Property


property_types = {
    "DATE_DRAW": "value_date",
    "DEMIS_ID_PC": "value_zip",
    "DEMIS_ID": "value_integer",
    "DESH_QC_PASSED": "value_varchar",
    "DESH_REJECTION_REASON": "value_varchar",
    "DOWNLOAD_ID": "value_varchar",
    "DOWNLOADING_TIMESTAMP": "value_varchar",
    "DUPLICATE_ID": "value_varchar",
    "HASHED_SEQUENCE": "value_varchar",
    "IMPORTED": "value_date",
    "LINEAGE": "value_varchar",
    "OWN_FASTA_ID": "value_varchar",
    "PROCESSING_DATE": "value_date",
    "PUBLICATION_STATUS": "value_varchar",
    "RECEIVE_DATE": "value_date",
    "SAMPLE_TYPE": "value_varchar",
    "SENDING_LAB_PC": "value_zip",
    "SENDING_LAB": "value_integer",
    "SEQ_REASON": "value_varchar",
    "SEQ_TYPE": "value_varchar",
    "STUDY": "value_varchar",
    "TIMESTAMP": "value_varchar",
    "VERSION": "value_integer",
}


def add_dummy_property(sample, property_name):
    property_value = get_dummy_value(property_name)
    property_type = property_types[property_name]
    property = Property.objects.get_or_create(name=property_name, datatype=property_type)[0]
    sample2property = Sample2Property.objects.create(
        sample=sample, property=property, **{property_type: property_value}
    )
    sample2property.save()


def get_dummy_value(property_name):
    type = property_types[property_name]
    if type == "value_varchar":
        return f"dummy_{property_name}_{str(random.randint(0, 1000000))}"
    elif type == "value_date":
        return f"18{str(random.randint(10, 99))}-{str(random.randint(1, 12)).zfill(2)}-{str(random.randint(1, 28)).zfill(2)}"
    elif type == "value_integer":
        return random.randint(0, 1000000)
    elif type == "value_zip":
        return f"dummy_{property_name}_{str(random.randint(0, 1000000))}.zip"


def fill_every_2nd_sample_dummy_properties(apps, schema_editor):
    samples = Sample.objects.all()
    if samples:
        for sample in samples[::2]:
            for property_name in property_types.keys():
                add_dummy_property(sample, property_name)


class Migration(migrations.Migration):
    dependencies = [
        ("rest_api", "0006_add_sample_dummy_properties"),
    ]

    operations = [
        migrations.RunPython(fill_every_2nd_sample_dummy_properties),
    ]
