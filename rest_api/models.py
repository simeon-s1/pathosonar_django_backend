from django.db import models


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


class Alignment2Variant(models.Model):
    alignment = models.ForeignKey(Alignment, models.DO_NOTHING)
    variant = models.ForeignKey("Variant", models.DO_NOTHING)

    class Meta:
        db_table = "alignment2variant"
        unique_together = (("variant", "alignment"),)


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
    standard = models.BigIntegerField(blank=True, null=True)
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
    type = models.CharField(blank=True, null=True)
    accession = models.CharField(unique=True, blank=True, null=True)
    symbol = models.CharField(blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    length = models.BigIntegerField(blank=True, null=True)
    segment = models.BigIntegerField(blank=True, null=True)
    standard = models.BigIntegerField(blank=True, null=True)

    class Meta:
        db_table = "molecule"


class Reference(models.Model):
    id = models.BigAutoField(primary_key=True)
    accession = models.CharField(unique=True, blank=True, null=True)
    description = models.CharField(blank=True, null=True)
    organism = models.CharField(blank=True, null=True)
    translation_group = models.ForeignKey("TranslationGroup", models.DO_NOTHING)
    standard = models.BigIntegerField(blank=True, null=True)

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
        max_digits=65535, decimal_places=65535, blank=True, null=True
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


class Variant(models.Model):
    id = models.BigAutoField(primary_key=True)
    element = models.ForeignKey(Element, models.DO_NOTHING, blank=True, null=True)
    ref = models.CharField(blank=True, null=True)
    alt = models.CharField(blank=True, null=True)
    start = models.BigIntegerField(blank=True, null=True)
    end = models.BigIntegerField(blank=True, null=True)
    parent_id = models.BigIntegerField(blank=True, null=True)
    label = models.CharField(blank=True, null=True)
    frameshift = models.BigIntegerField(blank=True, null=True)
    sample_name = models.CharField(max_length=50, blank=True, null=True)
    nuc_profile = models.CharField(max_length=1024, blank=True, null=True)
    aa_profile = models.CharField(max_length=1024, blank=True, null=True)
    imported = models.CharField(max_length=50, blank=True, null=True)
    collection_date = models.CharField(max_length=50, blank=True, null=True)
    release_date = models.CharField(max_length=50, blank=True, null=True)
    isolate = models.CharField(max_length=50, blank=True, null=True)
    length = models.IntegerField(blank=True, null=True)
    seq_tech = models.CharField(max_length=50, blank=True, null=True)
    country = models.CharField(max_length=50, blank=True, null=True)
    geo_location = models.CharField(max_length=50, blank=True, null=True)
    host = models.CharField(max_length=50, blank=True, null=True)
    genome_completeness = models.CharField(max_length=50, blank=True, null=True)
    reference_accession = models.CharField(max_length=50, blank=True, null=True)
    alignments = models.ManyToManyField(
        Alignment, through="Alignment2Variant", related_name="variants"
    )

    class Meta:
        db_table = "variant"


class Variant2Annotation(models.Model):
    variant = models.ForeignKey(Variant, models.DO_NOTHING)
    alignment = models.ForeignKey(Alignment, models.DO_NOTHING)
    annotation = models.ForeignKey(AnnotationType, models.DO_NOTHING)

    class Meta:
        db_table = "variant2annotation"
        unique_together = (("variant", "alignment", "annotation"),)


class Variant2Property(models.Model):
    property = models.ForeignKey(Property, models.DO_NOTHING)
    variant = models.ForeignKey(Variant, models.DO_NOTHING)
    value_integer = models.BigIntegerField(blank=True, null=True)
    value_float = models.DecimalField(
        max_digits=65535, decimal_places=65535, blank=True, null=True
    )
    value_text = models.TextField(blank=True, null=True)
    value_varchar = models.CharField(blank=True, null=True)
    value_blob = models.BinaryField(blank=True, null=True)
    value_date = models.DateField(blank=True, null=True)
    value_zip = models.CharField(blank=True, null=True)

    class Meta:
        db_table = "variant2property"
        unique_together = (
            ("property", "variant"),
            ("property", "variant"),
        )
