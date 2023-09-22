from rest_framework import serializers
from . import models

def find_or_create(data, model, serializer_class):
    try:
        return model.objects.get(**data)
    except model.DoesNotExist:
        serializer = serializer_class(data=data)
        serializer.is_valid(raise_exception=True)
        return serializer.save()

class MutationSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = models.Mutation
        fields = ["id", "label", "ref"]


class SamplePropertySerializer(serializers.ModelSerializer):
    value = serializers.SerializerMethodField()
    name = serializers.ReadOnlyField(source="property.name")

    class Meta:
        model = models.Sample2Property
        fields = ["value", "name"]

    def get_value(self, obj: models.Sample2Property):
        return (
            obj.value_integer
            or obj.value_float
            or obj.value_text
            or obj.value_varchar
            or obj.value_blob
            or obj.value_date
            or obj.value_zip
        )


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Sequence
        fields = "__all__"


class SampleSerializer(serializers.ModelSerializer):
    properties = SamplePropertySerializer(many=True, read_only=True)
    sequence = SequenceSerializer()

    def create(self, validated_data):
        validated_data["sequence"] = models.Sequence.objects.get_or_create(
            seqhash=validated_data["sequence"]
        )
        return super().create(validated_data)

    class Meta:
        model = models.Sample
        fields = "__all__"


class ReferenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Reference
        fields = "__all__"


class MoleculeSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Molecule
        fields = "__all__"


class ElementSerializer(serializers.ModelSerializer):
    def create(self, validated_data):
        validated_data["sequence"] = (
            validated_data["sequence"].strip().upper().replace("U", "T")
        )
        return super().create(validated_data)

    class Meta:
        model = models.Element
        fields = "__all__"


class ElempartSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Elempart
        fields = "__all__"
