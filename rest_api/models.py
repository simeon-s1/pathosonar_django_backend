from django.db import models
from django.db.models import UniqueConstraint


class Sequence(models.Model):
    id = models.BigAutoField(primary_key=True)
    seqhash = models.CharField(unique=True, max_length=200)

    class Meta:
        db_table = "sequence"


class Alignment(models.Model):
    id = models.BigAutoField(primary_key=True)
    element = models.BigIntegerField(blank=True, null=True)
    sequence = models.ForeignKey(
        "Sequence", models.DO_NOTHING, blank=True, null=True, related_name="alignments"
    )

    class Meta:
        db_table = "alignment"


class Alignment2Mutation(models.Model):
    alignment = models.ForeignKey(Alignment, models.DO_NOTHING)
    mutation = models.ForeignKey("Mutation", models.DO_NOTHING, blank=True, null=True)

    class Meta:
        db_table = "alignment2mutation"
        unique_together = (("mutation", "alignment"),)


class AnnotationType(models.Model):
    id = models.BigAutoField(primary_key=True)
    seq_ontology = models.CharField(max_length=50, blank=True, null=True)
    region = models.CharField(max_length=50, blank=True, null=True)
    test = models.CharField(max_length=50, blank=True, null=True)

    class Meta:
        db_table = "annotation_type"


class Element(models.Model):
    id = models.BigAutoField(primary_key=True)
    type = models.CharField(blank=True, null=True)
    accession = models.CharField(unique=True, blank=True, null=True)
    symbol = models.CharField(blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    start = models.BigIntegerField(blank=True, null=True)
    end = models.BigIntegerField(blank=True, null=True)
    strand = models.BigIntegerField(blank=True, null=True)
    sequence = models.TextField(blank=True, null=True)
    standard = models.BooleanField(blank=True, null=True)
    parent_id = models.BigIntegerField(blank=True, null=True)
    molecule = models.ForeignKey("Molecule", models.DO_NOTHING, blank=True, null=True)

    class Meta:
        db_table = "element"


class Elempart(models.Model):
    id = models.BigAutoField(primary_key=True)
    element = models.ForeignKey(Element, models.DO_NOTHING, blank=True, null=True)
    start = models.BigIntegerField()
    end = models.BigIntegerField()
    strand = models.BigIntegerField()
    base = models.FloatField()
    segment = models.BigIntegerField()

    class Meta:
        db_table = "elempart"


class Lineages(models.Model):
    id = models.BigAutoField(primary_key=True)
    lineage = models.CharField(max_length=100)
    sublineage = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "lineages"


class Molecule(models.Model):
    id = models.BigAutoField(primary_key=True)
    reference = models.ForeignKey("Reference", models.DO_NOTHING, blank=True, null=True)
    accession = models.CharField(unique=True, blank=True, null=True)
    symbol = models.CharField(blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    length = models.BigIntegerField(blank=True, null=True)
    segment = models.BigIntegerField(blank=True, null=True)
    standard = models.BooleanField(blank=True, null=True)
    type = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "molecule"


class Reference(models.Model):
    id = models.BigAutoField(primary_key=True)
    accession = models.CharField(unique=True, blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    organism = models.CharField(blank=True, null=True)
    translation_group = models.ForeignKey("TranslationGroup", models.DO_NOTHING)
    standard = models.BooleanField(blank=True, null=True)

    class Meta:
        db_table = "reference"


class Property(models.Model):
    id = models.BigAutoField(primary_key=True)
    name = models.CharField(unique=True, blank=True, null=True)
    datatype = models.CharField(blank=True, null=True)
    querytype = models.CharField(blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    target = models.CharField(blank=True, null=True)
    standard = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "property"


class Sample(models.Model):
    id = models.BigAutoField(primary_key=True)
    name = models.CharField(unique=True, blank=True, null=True)
    datahash = models.CharField(blank=True, null=True)
    sequence = models.ForeignKey(
        Sequence, models.DO_NOTHING, blank=True, null=True, related_name="samples"
    )

    class Meta:
        db_table = "sample"


class Sample2Property(models.Model):
    property = models.ForeignKey(Property, models.DO_NOTHING)
    sample = models.ForeignKey(Sample, models.DO_NOTHING, related_name="properties")
    value_integer = models.BigIntegerField(blank=True, null=True)
    value_float = models.DecimalField(
        max_digits=10, decimal_places=10, blank=True, null=True
    )
    value_text = models.TextField(blank=True, null=True)
    value_varchar = models.CharField(blank=True, null=True)
    value_blob = models.BinaryField(blank=True, null=True)
    value_date = models.DateField(blank=True, null=True)
    value_zip = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "sample2property"
        unique_together = (("property", "sample"),)


class Translation(models.Model):
    id = models.BigAutoField(primary_key=True)
    group = models.ForeignKey(
        "TranslationGroup", models.DO_NOTHING, blank=True, null=True
    )
    codon = models.CharField(blank=True, null=True)
    aa = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "translation"
        unique_together = (("group", "codon"),)


class TranslationGroup(models.Model):
    id = models.BigAutoField(primary_key=True)

    class Meta:
        db_table = "translation_group"


class Mutation(models.Model):
    id = models.BigAutoField(primary_key=True)
    element = models.ForeignKey(Element, models.DO_NOTHING, blank=True, null=True)
    ref = models.CharField(blank=True, null=True)
    alt = models.CharField(blank=True, null=True)
    start = models.BigIntegerField(blank=True, null=True)
    end = models.BigIntegerField(blank=True, null=True)
    parent_id = models.BigIntegerField(blank=True, null=True)
    label = models.CharField(blank=True, null=True)
    frameshift = models.BigIntegerField(blank=True, null=True)
    alignments = models.ManyToManyField(
        Alignment, through="Alignment2Mutation", related_name="mutations"
    )

    class Meta:
        db_table = "mutation"
        UniqueConstraint(
            name="unique_mutation",
            fields=[
                "element",
                "ref",
                "alt",
                "start",
                "end",
                "parent_id",
                "label",
                "frameshift",
            ]
        )


class Mutation2Annotation(models.Model):
    mutation = models.ForeignKey(Mutation, models.DO_NOTHING, blank=True, null=True)
    alignment = models.ForeignKey(Alignment, models.DO_NOTHING, blank=True, null=True)
    annotation = models.ForeignKey(
        AnnotationType, models.DO_NOTHING, blank=True, null=True
    )

    class Meta:
        db_table = "mutation2annotation"
        unique_together = (("mutation", "alignment", "annotation"),)


class Mutation2Property(models.Model):
    property = models.ForeignKey(Property, models.DO_NOTHING)
    mutation = models.ForeignKey(Mutation, models.DO_NOTHING)
    value_integer = models.BigIntegerField(blank=True, null=True)
    value_float = models.DecimalField(
        max_digits=10, decimal_places=10, blank=True, null=True
    )
    value_text = models.TextField(blank=True, null=True)
    value_varchar = models.CharField(blank=True, null=True)
    value_blob = models.BinaryField(blank=True, null=True)
    value_date = models.DateField(blank=True, null=True)
    value_zip = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "mutation2property"
        unique_together = (
            ("property", "mutation"),
            ("property", "mutation"),
        )


class EnteredData(models.Model):
    type = models.CharField(max_length=50, blank=True, null=True)
    name = models.CharField(max_length=400, blank=True, null=True)
    date = models.DateField(blank=True, null=True)

    class Meta:
        db_table = "entered_data"
        unique_together = (("type", "name"),)
