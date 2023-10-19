import pathlib
import pickle
from dataclasses import dataclass
from rest_api.models import Sequence, Sample, Alignment, Mutation, Gene
from rest_api.serializers import SampleSerializer, MutationSerializer, find_or_create
from threading import Lock


@dataclass
class SampleRaw:
    algn_file: str
    algnid: int | None
    anno_tsv_file: str
    anno_vcf_file: str
    cds_file: str
    header: str
    lift_file: str
    mafft_seqfile: str
    name: str
    properties: dict
    ref_file: str
    refmol: str
    refmolid: int
    sampleid: int | None
    seq_file: str
    seqhash: str
    sourceid: int
    translationid: int
    tt_file: str
    var_file: str
    vcffile: str


@dataclass
class VarRaw:
    ref: str
    start: int
    end: int
    alt: str | None
    element_id: int
    label: str
    frameshift: bool


class SampleImport:
    def __init__(
        self,
        path: pathlib.Path,
        mutation_cache: dict[str, Mutation],
        import_folder="import_data",
    ):        
        self.sample_file_path = path
        self.import_folder = import_folder
        self.mutation_cache = mutation_cache
        self.sample_raw = SampleRaw(**self._import_pickle(path))
        self.vars_raw = [var for var in self._import_vars(self.sample_raw.var_file)]
        self.seq = "".join(
            [line for line in self._import_seq(self.sample_raw.seq_file)]
        )

    def write_to_db(self):
        sequence = Sequence.objects.get_or_create(seqhash=self.sample_raw.seqhash)[0]
        sample = self._find_or_create_sample(sequence)
        gene = Gene.objects.get(id=1)
        alignment = Alignment.objects.get_or_create(sequence=sequence, gene=gene)[0]
        for var_raw in self.vars_raw:
            if self.mutation_cache.get(var_raw.label):
                mutation = self.mutation_cache[var_raw.label]
            else:
                gene = Gene.objects.get(id=2)
                mutation_data = {
                    "gene": gene.id,
                    "ref": var_raw.ref,
                    "alt": var_raw.alt,
                    "start": var_raw.start,
                    "end": var_raw.end,
                    "label": var_raw.label,
                    "frameshift": int(var_raw.frameshift),
                }
                mutation = find_or_create(mutation_data, Mutation, MutationSerializer)
                self.mutation_cache[var_raw.label] = mutation
            mutation.alignments.add(alignment)
            mutation.save()
        return sample

    def _find_or_create_sample(self, sequence):
        try:
            return Sample.objects.get(name=self.sample_raw.name, sequence=sequence)
        except Sample.DoesNotExist:
            sample_serializer = SampleSerializer(
                data={
                    "name": self.sample_raw.name,
                    "sequence": sequence.id,
                    "properties": self.sample_raw.properties,
                }
            )
            sample_serializer.is_valid(raise_exception=True)
            return sample_serializer.save()

    def _import_pickle(self, path: str):
        with open(path, "rb") as f:
            return pickle.load(f)

    def _import_vars(self, path):
        file_name = pathlib.Path(path).name
        self.var_file_path = (
            pathlib.Path(self.import_folder)
            .joinpath("var")
            .joinpath(file_name[:2])
            .joinpath(file_name)
        )
        with open(self.var_file_path, "r") as handle:
            for line in handle:
                if line == "//":
                    break
                var_raw = line.strip("\r\n").split("\t")
                yield VarRaw(
                    var_raw[0],  # ref
                    int(var_raw[1]),  # start
                    int(var_raw[2]),  # end
                    None if var_raw[3] == " " else var_raw[3],  # alt
                    int(var_raw[4]),  # element id
                    var_raw[5],  # label
                    bool(var_raw[6]),  # frameshift
                )

    def _import_seq(self, path):
        file_name = pathlib.Path(path).name
        self.seq_file_path = (
            pathlib.Path(self.import_folder)
            .joinpath("seq")
            .joinpath(file_name[:2])
            .joinpath(file_name)
        )
        with open(self.seq_file_path, "r") as handle:
            for line in handle:
                yield line.strip("\r\n")

    def move_files(self, success):
        for file in [
            self.sample_file_path,
            self.var_file_path,
            self.seq_file_path,
        ]:
            import_data_location = pathlib.Path(file).parent.parent.parent
            folder_structure = pathlib.Path(file).relative_to(import_data_location)
            if success:
                target_folder = pathlib.Path(import_data_location).joinpath(
                    "imported_successfully"
                )
            else:
                target_folder = pathlib.Path(import_data_location).joinpath(
                    "failed_import"
                )
            target_folder.joinpath(folder_structure.parent).mkdir(
                parents=True, exist_ok=True
            )
            pathlib.Path(file).rename(target_folder.joinpath(folder_structure))