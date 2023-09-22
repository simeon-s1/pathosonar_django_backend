from Bio import SeqFeature, SeqIO, SeqRecord
from django.core.files.uploadedfile import InMemoryUploadedFile

from rest_api.models import Element, Elempart, Molecule, Reference
from rest_api.serializers import (
    ElementSerializer,
    ElempartSerializer,
    MoleculeSerializer,
    ReferenceSerializer,
    find_or_create,
)


def import_gbk_file(uploaded_file: InMemoryUploadedFile, standard=False):
    records = list(SeqIO.parse(_temp_save_file(uploaded_file), "genbank"))
    records: list[SeqRecord.SeqRecord]
    reference = _put_reference_from_record(records[0], standard)
    standard = 1
    for i, record in enumerate(records):
        molecule_data = {
            "reference": reference.id,
            "accession": f"{record.name}.{record.annotations['sequence_version']}",
            "symbol": record.annotations.get("symbol", ""),
            "description": record.description,
            "segment": i,
            "standard": standard,  # TODO
        }
        molecule = find_or_create(molecule_data, Molecule, MoleculeSerializer)

        source_features = list(filter(lambda x: x.type == "source", record.features))
        if len(source_features) != 1:
            raise ValueError("Expecting exactly one source feature.")
        source_feature = source_features[0]
        source_feature: SeqFeature.SeqFeature
        if source_feature.location is None:
            raise ValueError("No location information found for source feature.")
        source_element_data = {
            "molecule": molecule.id,
            "type": "source",
            "start": source_feature.location.start.position,
            "end": source_feature.location.end.position,
            "sequence": str(source_feature.extract(record.seq)),
            "standard": "1",
        }
        molecule.type = source_feature.qualifiers.get("mol_type", [""])[0]
        molecule.length = len(source_element_data["sequence"])
        # molecule.segment = source_feature.qualifiers.get("segment", [""])[0]
        molecule.save()

        source_element = find_or_create(source_element_data, Element, ElementSerializer)
        _create_elemparts(source_feature, source_element)
        for feature in record.features:
            if "pseudogene" in feature.qualifiers or feature.type not in [
                "CDS",
                "gene",
            ]:
                continue
            element = _put_element_from_feature(feature, record.seq, source_element.id)
            _create_elemparts(feature, element)
        standard = 0
    return records


def _process_segments(
    feat_location_parts: list[SeqFeature.FeatureLocation | SeqFeature.CompoundLocation],
    cds: bool = False,
) -> list[dict[str, int]]:
    """
    Process the genomic regions (segments) of a feature.

    Args:
        feat_location_parts (List[Union[FeatureLocation, CompoundLocation]]): List of feature location parts.
        cds (bool): A flag indicating whether the segment corresponds to a coding sequence.
                    Default is False.

    Returns:
        segments (List[List[int]]): A list of processed segments. Each segment is represented
                                    as a list of integers [start, end, strand, base, index].
    """
    base = 0
    div = 1 if not cds else 3
    segments = []
    for i, segment in enumerate(feat_location_parts, 1):
        segments.append(
            {
                "start": segment.start.position,
                "end": segment.end.position,
                "strand": segment.strand,
                "base": base,
                "segment": i,
            }
        )
        base += round((segment.end.position - segment.start.position - 1) / div, 1)
    return segments


def _validate_segment_lengths(parts, accession):
    parts = _process_segments(parts, cds=True)
    if sum([abs(x["end"] - x["start"]) for x in parts]) % 3 != 0:
        raise ValueError(f"The length of cds '{accession}' is not a multiple of 3.")


def _determine_accession(feature: SeqFeature.SeqFeature) -> str | None:
    if feature.id != "<unknown id>":
        return feature.id
    elif "gene" in feature.qualifiers and feature.type == "gene":
        return feature.qualifiers["gene"][0]
    elif "protein_id" in feature.qualifiers and feature.type == "CDS":
        return feature.qualifiers["protein_id"][0]
    elif "locus_tag" in feature.qualifiers:
        return feature.qualifiers["locus_tag"][0]


def _put_element_from_feature(
    feature: SeqFeature.SeqFeature, source_seq: str, source_id: int
) -> Element:
    if feature.type not in ["gene", "CDS"]:
        raise ValueError(f"The provided feature type({feature.type}) is unknown.")
    if feature.location is None:
        raise ValueError("No location information found for gene feature.")

    element_data = {
        "type": feature.type.lower(),
        "start": feature.location.start.position,
        "end": feature.location.end.position,
        "strand": feature.strand,
        "parent_id": source_id,
        "standard": 0,
    }

    if accession := _determine_accession(feature):
        element_data["accession"] = accession
    else:
        raise ValueError("No qualifier for gene accession found.")

    if "gene" in feature.qualifiers:
        element_data["symbol"] = feature.qualifiers["gene"][0]
    elif "locus_tag" in feature.qualifiers:
        element_data["symbol"] = feature.qualifiers["locus_tag"][0]
    else:
        raise ValueError("No qualifier for gene symbol found.")

    if feature.type == "CDS":
        _validate_segment_lengths(feature.location.parts, element_data["accession"])
        element_data["sequence"] = feature.qualifiers.get("translation", [""])[0]
        element_data["description"] = feature.qualifiers.get("product", [""])[0]
    else:
        element_data["sequence"] = str(feature.extract(source_seq))
        element_data["description"] = ""
    return find_or_create(element_data, Element, ElementSerializer)


def _put_reference_from_record(
    record: SeqRecord.SeqRecord, standard: bool
) -> Reference:
    reference = {
        "accession": f"{record.name}.{record.annotations['sequence_version']}",
        "description": record.description,
        "organism": record.annotations["organism"],
        "translation_group": 1,
        "standard": standard,
    }
    return find_or_create(reference, Reference, ReferenceSerializer)


def _create_elemparts(feature: SeqFeature.SeqFeature, element: Element):
    for elempart in _process_segments(feature.location.parts):
        elempart_data = {
            "element": element.id,
            **elempart,
        }
        find_or_create(elempart_data, Elempart, ElempartSerializer)


def _temp_save_file(uploaded_file: InMemoryUploadedFile):
    file_path = pathlib.Path("import_data") / uploaded_file.name
    with open(file_path, "wb") as f:
        f.write(uploaded_file.read())
    return file_path
