from rest_framework import serializers
from . import models


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