import pathlib
import pickle
from datetime import datetime

from Bio import SeqIO, SeqRecord, SeqFeature

from rest_api.models import (
    EnteredData,
    Mutation,
    Reference,
    Molecule,
    Element,
    Elempart,
)
from rest_api.serializers import (
    find_or_create,
    MutationSerializer,
    ReferenceSerializer,
    SampleSerializer,
    MoleculeSerializer,
    ElementSerializer,
    ElempartSerializer,
)
from django.core.files.uploadedfile import InMemoryUploadedFile

# tsv format

mutation_tsv_header = [
    "ref",
    "start",
    "end",
    "alt",
    "parent_id",
    "label",
    "frameshift",
]

extension_to_model = {".var": Mutation}

extension_to_header_map = {".var": mutation_tsv_header}

extension_to_serializer_map = {".var": MutationSerializer, ".sample": SampleSerializer}

pickle_extensions = [".sample"]

import_order = [".var", ".sample"]


def run_data_entry():
    print("run_data_entry")
    for extension in import_order:
        files = pathlib.Path("import_data").glob("**/*")
        objects = []
        for file in files:
            if file.is_file() and file.suffix == extension:
                base_name = file.name
                if not EnteredData.objects.filter(
                    type=extension, name=base_name
                ).exists():
                    try:
                        print("import ", file)
                        import_file(file, extension)
                    except Exception as e:
                        print("Error importing file: ", base_name, e)
                        raise e
                    else:
                        EnteredData.objects.create(
                            type=extension, name=base_name, date=datetime.now()
                        )
                    finally:
                        print("import done")
                else:
                    print("already imported ", file)
        serializer_class = extension_to_serializer_map[extension]
        serializer = serializer_class(data=objects, many=True)
        serializer.is_valid(raise_exception=True)
        serializer.save()


def import_file(path, extension):
    if extension in pickle_extensions:
        objects = import_pickle(path, extension)
    else:
        objects = import_tsv(path, extension)


def import_pickle(path, extension):
    with open(path, "rb") as f:
        objects = pickle.load(f)
        print("objects", objects)
    raise Exception("not implemented")
    return objects


def import_tsv(path, extension):
    objects = []
    with open(path) as f:
        if extension not in extension_to_header_map:
            raise Exception("Unknown file extension: " + extension)
        header = extension_to_header_map[extension]
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith(r"//"):
                continue
            fields = line.split("\t")
            if len(fields) != len(header):
                raise Exception(
                    "Expected "
                    + str(len(header))
                    + " fields, got "
                    + str(len(fields))
                    + " fields in line: "
                    + line
                )
            d = dict(zip(header, fields))
            objects.append(d)
    return objects

