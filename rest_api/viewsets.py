from dataclasses import dataclass
import json
import time
from django.db.models import F, Count, QuerySet, Prefetch
from rest_framework import serializers, viewsets, generics
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.request import Request
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rest_api.models.sample import Sample
    from rest_api.models.variant import Variant

from . import models


class ElementSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = models.Element
        fields = [
            "id",
            "accession",
            "symbol",
            "description",
            "start",
            "end",
            "strand",
            "sequence",
            "standard",
            "parent",
        ]


class VariantSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = models.Variant
        fields = ["id", "label", "ref"]


class ElementReferencesSerializer(serializers.HyperlinkedModelSerializer):
    id = serializers.ReadOnlyField(source="molecule.reference.id")
    accession = serializers.ReadOnlyField(source="molecule.reference.accession")

    class Meta:
        model = models.Element
        fields = ["id", "accession", "sequence"]


class ElementViewSet(viewsets.ModelViewSet):
    queryset = models.Element.objects.all()
    serializer_class = ElementSerializer

    @action(detail=False, methods=["get"])
    def references(self, request: Request, *args, **kwargs):
        queryset = models.Element.objects.filter(type="source")
        serializers = ElementReferencesSerializer(queryset, many=True)
        return Response(serializers.data)

    @action(detail=False, methods=["get"])
    def distinct_genes(self, request: Request, *args, **kwargs):
        queryset = models.Element.objects.distinct("symbol").values("symbol")
        if ref := request.query_params.get("reference"):
            queryset = queryset.filter(molecule__reference__accession=ref)
        return Response({"genes": [item["symbol"] for item in queryset]})


class VariantViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Variant.objects.all()
    serializer_class = VariantSerializer


class PropertySerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Property
        fields = "__all__"
        filter_backends = [DjangoFilterBackend]


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


class SampleSerializer(serializers.ModelSerializer):
    properties = SamplePropertySerializer(many=True, read_only=True)

    class Meta:
        model = models.Sample
        fields = "__all__"


class SampleViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Sample.objects.all().order_by("id")
    serializer_class = SampleSerializer
    filter_backends = [DjangoFilterBackend]

    @action(detail=False, methods=["get"])
    def count_unique_nt_mut_ref_view_set(self, request: Request, *args, **kwargs):
        # TODO-smc abklären ob das so richtig ist
        queryset = (
            models.Variant.objects.exclude(
                element__type="cds",
                alt="N",
            )
            .values("element__molecule__reference__accession")
            .annotate(count=Count("id"))
        )
        dict = {
            item["element__molecule__reference__accession"]: item["count"]
            for item in queryset
        }
        return Response(data=dict)

    @action(detail=False, methods=["get"])
    def genomes(self, request: Request, *args, **kwargs):
        queryset = models.Sample.objects.all()
        if property_filters := request.query_params.get("properties"):
            property_filters = json.loads(property_filters)
            for property_name, filter in property_filters.items():
                datatype = models.Property.objects.get(name=property_name).datatype
                query = {f"properties__property__name": property_name}
                for key, value in filter.items():
                    query[f"properties__{datatype}__{key}"] = value
                queryset = queryset.filter(**query)
        profile_filter_methods = {
            "snpProfileNtFilters": self.filter_snp_profile_nt,
            "snpProfileAAFilters": self.filter_snp_profile_aa,
            "delProfileNtFilters": self.filter_del_profile_nt,
            "delProfileAAFilters": self.filter_del_profile_aa,
            "insProfileNtFilters": self.filter_ins_profile_nt,
            "insProfileAAFilters": self.filter_ins_profile_aa,
        }
        for (
            profile_filter_name,
            profile_filter_method,
        ) in profile_filter_methods.items():
            if profile_filters := request.query_params.get(profile_filter_name):
                profile_filters = json.loads(profile_filters)
                for profile_filter in profile_filters:
                    queryset = profile_filter_method(**profile_filter, qs=queryset)
        if snp_profile_filters := request.query_params.get("snp_profile_filters"):
            snp_profile_filters = json.loads(snp_profile_filters)
            for snp_profile_filter in snp_profile_filters:
                queryset = self.filter_snp_profile_nt(**snp_profile_filter, qs=queryset)
        queryset = queryset.prefetch_related("properties__property")
        # obj.sequence.alignments.filter(variants__frameshift=1).values_list("variants__label", flat=True)
        queryset = queryset.prefetch_related(
            Prefetch(
                "sequence__alignments__variants",
                queryset=models.Variant.objects.filter(frameshift="1"),
                to_attr="frameshifts",
            )
        )

        queryset = self.paginate_queryset(queryset)
        serializer = SampleGenomesSerializer(queryset, many=True)
        data = serializer.data
        return self.get_paginated_response(data)

    @action(detail=False, methods=["get"])
    def download_genomes_export(self, request: Request, *args, **kwargs):
        start_time = time.time()
        queryset = models.Sample.objects.all()
        if property_filters := request.query_params.get("properties"):
            property_filters = json.loads(property_filters)
            for property_name, filter in property_filters.items():
                datatype = models.Property.objects.get(name=property_name).datatype
                query = {f"properties__property__name": property_name}
                for key, value in filter.items():
                    query[f"properties__{datatype}__{key}"] = value
                queryset = queryset.filter(**query)
        queryset.prefetch_related(
            "sequence__alignments__variants__element__molecule__reference"
        )
        queryset.prefetch_related("properties__property")
        # save queryset as file:
        line_count = 0
        with open("genomes_export.csv", "w") as file:
            for sample in queryset.iterator():
                file.writelines(str(SampleGenomesSerializer(sample).data))
                file.write("\n")
                line_count += 1
                if line_count % 100 == 0:
                    print(f"{line_count} lines written")
        print("--- %s seconds ---" % (time.time() - start_time))
        return Response("OK")

    def filter_snp_profile_nt(
        self,
        ref_nuc: str,
        ref_pos: int,
        alt_nuc: str,
        qs: QuerySet | None = None,
        exclude: bool = False,
    ) -> QuerySet:
        # For NT: ref_nuc followed by ref_pos followed by alt_nuc (e.g. T28175C).
        if qs is None:
            qs = Variant.objects.all()
        filters = {
            "start": ref_pos,
            "ref": ref_nuc,
            "alt": alt_nuc,
        }
        qs = qs.exclude(**filters) if exclude else qs.filter(**filters)
        return qs

    def filter_snp_profile_aa(
        self,
        protein_symbol: str,
        ref_aa: str,
        ref_pos: str,
        alt_aa: str,
        qs: QuerySet | None = None,
        exclude: bool = False,
    ) -> QuerySet:
        # For AA: protein_symbol:ref_aa followed by ref_pos followed by alt_aa (e.g. OPG098:E162K)
        if qs is None:
            qs = Variant.objects.filter(element__molecule__symbol=protein_symbol)
        filters = {
            "start": ref_pos,
            "ref": ref_aa,
            "alt": alt_aa,
        }
        qs = qs.exclude(**filters) if exclude else qs.filter(**filters)
        return qs

    def filter_del_profile_nt(
        self,
        first_deleted_nt: str,
        last_deleted_nt: str,
        qs: QuerySet | None = None,
    ) -> QuerySet:
        # For NT: del:first_NT_deleted-last_NT_deleted (e.g. del:133177-133186).
        if qs is None:
            qs = Variant.objects.all()

    def filter_del_profile_aa(
        self,
        protein_symbol: str,
        first_deleted_aa: int,
        last_deleted_aa: int,
        exclude: bool = False,
        qs: QuerySet | None = None,
    ) -> QuerySet:
        # For AA: protein_symbol:del:first_AA_deleted-last_AA_deleted (e.g. OPG197:del:34-35)
        if qs is None:
            qs = Variant.objects.filter(element__molecule__symbol=protein_symbol)
        filters = {"start": first_deleted_aa - 1, "end": last_deleted_aa, "alt": ""}
        qs = qs.exclude(**filters) if exclude else qs.filter(**filters)
        return qs

    def filter_ins_profile_nt(
        self,
        ref_nuc: str,
        ref_pos: int,
        alt_nucs: str,
        qs: QuerySet | None = None,
        exclude: bool = False,
    ) -> QuerySet:
        # For NT: ref_nuc followed by ref_pos followed by alt_nucs (e.g. T133102TTT)
        if qs is None:
            qs = Variant.objects.all()
        filters = {
            "start": ref_pos,
            "ref": ref_nuc,
            "alt": alt_nucs,
        }
        qs = qs.exclude(**filters) if exclude else qs.filter(**filters)
        return qs

    def filter_ins_profile_aa(
        self,
        protein_symbol: str,
        ref_aa: str,
        ref_pos: int,
        alt_aas: str,
        qs: QuerySet | None = None,
        exclude: bool = False,
    ) -> QuerySet:
        # For AA: protein_symbol:ref_aa followed by ref_pos followed by alt_aas (e.g. OPG197:A34AK)
        if qs is None:
            qs = Variant.objects.filter(element__molecule__symbol=protein_symbol)
        filters = {
            "start": ref_pos,
            "ref": ref_aa,
            "alt": alt_aas,
        }
        qs = qs.exclude(**filters) if exclude else qs.filter(**filters)
        return qs


class ReferenceSerializer(serializers.HyperlinkedModelSerializer):
    sequence = serializers.CharField(source="molecules.elements.sequence")

    class Meta:
        model = models.Reference
        fields = ["accession", "description", "organism", "standard", "sequence"]


class ReferenceViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Reference.objects.all()
    serializer_class = ReferenceSerializer


class GenesSerializer(serializers.HyperlinkedModelSerializer):
    reference_accession = serializers.CharField(source="molecule.reference.accession")

    class Meta:
        model = models.Element
        fields = [
            "reference_accession",
            "type",
            "symbol",
            "description",
            "start",
            "end",
            "strand",
            "sequence",
        ]


class GenesViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Element.objects.all()
    serializer_class = GenesSerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ["molecule__reference__accession", "type"]


class SNP1Serializer(serializers.HyperlinkedModelSerializer):
    reference_accession = serializers.CharField(
        source="element.molecule.reference.accession"
    )

    class Meta:
        model = models.Variant
        fields = [
            "reference_accession",
            "ref",
            "alt",
            "start",
            "end",
        ]


class SNP1ViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Variant.objects.filter(
        ref__in=["C", "T", "G", "A"], alt__in=["C", "T", "G", "A"]
    ).exclude(ref=F("alt"))
    serializer_class = SNP1Serializer


class MutationSignatureSerializer(serializers.HyperlinkedModelSerializer):
    reference_accession = serializers.CharField(
        source="element.molecule.reference.accession"
    )
    count = serializers.IntegerField()

    class Meta:
        model = models.Variant
        fields = [
            "reference_accession",
            "ref",
            "alt",
            "start",
            "end",
            "count",
        ]


class MutationSignatureRawSerializer(serializers.HyperlinkedModelSerializer):
    reference_accession = serializers.ReadOnlyField(
        source="element.molecule.reference.accession"
    )

    class Meta:
        model = models.Variant
        fields = [
            "reference_accession",
            "ref",
            "alt",
            "start",
            "end",
        ]


class MutationSignatureViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = (
        models.Variant.objects.filter(ref__in=["C", "T"], alt__in=["C", "T", "G", "A"])
        .exclude(ref=F("alt"))
        .annotate(count=Count("alignments__sequence__samples"))
    )
    serializer_class = MutationSignatureSerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ["element__molecule__reference__accession", "ref", "alt"]

    @action(detail=False, methods=["get"])
    def raw(self, request: Request, *args, **kwargs):
        queryset = models.Variant.objects.filter(
            ref__in=["C", "T"], alt__in=["C", "T", "G", "A"]
        ).exclude(ref=F("alt"))
        page = self.paginate_queryset(queryset)
        serializer = MutationSignatureRawSerializer(page, many=True)
        if page is not None:
            return self.get_paginated_response(serializer.data)
        return self.get_paginated_response(serializer.data)


class PropertyViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
    generics.mixins.RetrieveModelMixin,
):
    queryset = models.Sample2Property.objects.all()
    serializer_class = PropertySerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ["sample__id", "property__name"]

    @action(detail=False, methods=["get"])
    def unique_collection_dates(self, request: Request, *args, **kwargs):
        queryset = models.Sample2Property.objects.filter(
            property__name="COLLECTION_DATE", value_date__isnull=False
        ).distinct("value_date")
        if ref := request.query_params.get("reference"):
            queryset = queryset.filter(
                sample__sequence__alignment__element__molecule__reference__accession=ref
            )
        queryset = queryset.distinct("value_text")
        date_list = [item.value_date for item in queryset]
        return Response(data={"collection_dates": date_list})

    @action(detail=False, methods=["get"])
    def unique_countries(self, request: Request, *args, **kwargs):
        queryset = models.Sample2Property.objects.filter(
            property__name="COUNTRY", value_text__isnull=False
        )
        if ref := request.query_params.get("reference"):
            queryset = queryset.filter(
                sample__sequence__alignment__element__molecule__reference__accession=ref
            )
        queryset = queryset.distinct("value_text")
        country_list = [item.value_text for item in queryset]
        return Response(data={"countries": country_list})

    @action(detail=False, methods=["get"])
    def unique_sequencing_techs(self, request: Request, *args, **kwargs):
        queryset = models.Sample2Property.objects.filter(
            property__name="SEQ_TECH", value_text__isnull=False
        )
        if ref := request.query_params.get("reference"):
            queryset = queryset.filter(
                sample__sequence__alignment__element__molecule__reference__accession=ref
            )
        queryset = queryset.distinct("value_text")
        sequencing_tech_list = [item.value_text for item in queryset]
        return Response(data={"sequencing_techs": sequencing_tech_list})

    @action(detail=False, methods=["get"])
    def distinct_property_names(self, request: Request, *args, **kwargs):
        queryset = models.Property.objects.all()
        queryset = queryset.distinct("name")
        property_names = [item.name for item in queryset]
        return Response(data={"property_names": property_names})


class MutationFrequencySerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = models.Variant
        fields = [
            "element_symbol",
            "variant_label",
            "count",
        ]


class AAMutationViewSet(
    viewsets.GenericViewSet,
    generics.mixins.ListModelMixin,
):
    # input sequencing_tech list, country list, gene list, include partial bool,
    # reference_value int, min_nb_freq int = 1?
    queryset = models.Element.objects.all()

    @action(detail=False, methods=["get"])
    def mutation_frequency(self, request: Request, *args, **kwargs):
        country_list = request.query_params.getlist("countries")
        sequencing_tech_list = request.query_params.getlist("seq_techs")
        gene_list = request.query_params.getlist("genes")
        include_partial = bool(request.query_params.get("include_partial"))
        reference_value = request.query_params.get("reference_value")
        min_nb_freq = request.query_params.get("min_nb_freq")

        samples_query = models.Sample.objects.filter(
            sample2property__property__name="COUNTRY",
            sample2property__value_text__in=country_list,
        ).filter(
            sample2property__property__name="SEQ_TECH",
            sample2property__value_text__in=sequencing_tech_list,
        )
        if not include_partial:
            samples_query.filter(
                sample2property__property__name="GENOME_COMPLETENESS",
                sample2property__value_text="complete",
            )
        variant_query = (
            models.Variant.objects.filter(
                element__molecule__reference__accession=reference_value
            )
            .filter(alignments__sequence__samples__in=samples_query)
            .filter(element__symbol__in=gene_list)
            .annotate(variant_count=Count("alignments__sequence__samples"))
            .filter(variant_count__gte=min_nb_freq)
            .order_by("-variant_count")
        )
        response = [
            {
                "symbol": variant.element.symbol,
                "variant": variant.label,
                "count": variant.variant_count,
            }
            for variant in variant_query
        ]

        return Response(data=response)


class VariantLabelSerializer(serializers.Serializer):
    labels = serializers.ListField()


class SampleGenomesSerializer(serializers.HyperlinkedModelSerializer):
    properties = SamplePropertySerializer(many=True, read_only=True)
    frameshifts = serializers.SerializerMethodField()
    # genomic_profiles = serializers.SerializerMethodField()
    # proteomic_profiles = serializers.SerializerMethodField()

    class Meta:
        model = models.Sample
        fields = [
            "id",
            "name",
            "sequence_id",
            "datahash",
            "properties",
            "frameshifts",
            # "genomic_profiles",
            # "proteomic_profiles",
        ]

    def get_frameshifts(self, obj: models.Sample):
        return obj.frameshifts

    def get_genomic_profiles(self, obj: models.Sample):
        query = (
            obj.sequence.alignments.exclude(variants__alt="N")
            .filter(variants__element__type="gene")
            .values_list("variants__label", flat=True)
        )
        return query

    def get_proteomic_profiles(self, obj: models.Sample):
        query = (
            obj.sequence.alignments.exclude(variants__alt="X")
            .filter(variants__element__type="cds")
            .values_list("variants__label", flat=True)
        )
        return VariantLabelSerializer({"labels": query}).data

    def csv_line(self):
        line = []
        for key in self.data.keys():
            if key == "properties":
                for prop in self.data[key]:
                    line.append(str(prop["value"]))
            else:
                line.append(str(self.data[key]))
        return ",".join(line)


class SampleGenomeViewSet(viewsets.GenericViewSet, generics.mixins.ListModelMixin):
    queryset = models.Sample.objects.all()
    serializer_class = SampleGenomesSerializer

    @action(detail=False, methods=["get"])
    def match(self, request: Request, *args, **kwargs):
        profile_filters = request.query_params.getlist("profile_filters")
        param_filters = request.query_params.getlist("param_filters")

    @action(detail=False, methods=["get"])
    def test_profile_filters():
        pass
