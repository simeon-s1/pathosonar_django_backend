import pathlib
from datetime import datetime
from multiprocessing import Lock, Pool, cpu_count
import traceback

import django
from django.db import connections, transaction

from rest_api.data_entry.sample_import import SampleImport
from rest_api.models import EnteredData

SAMPLE_EXTENSION = ".sample"


log_lock = Lock()
debug = False  # set to True to run only one file
_async = False


def async_init():
    django.setup()
    connections.close_all()


class SampleEntryJob:
    def __init__(self):
        self.start_timer = datetime.now()
        self.files_number = 0
        self.files = []
        self.error_log_file = open(
            f"{self.start_timer.strftime('%Y%m%d-%H%M%S')}_sample_import_error_log.txt",
            "w",
        )
        self.return_list = []
        self.list_lock = Lock()

    def run_data_entry(self):
        print("--- running data entry ---")
        self.files = list(
            pathlib.Path("import_data").joinpath("samples").glob("**/*.sample")
        )
        self.files_number = len(self.files)
        self.last_estimate_printed_at = 0
        print(f"{self.files_number} files found")
        timer = datetime.now()
        if debug:
            for file in self.files:
                if file.suffix == SAMPLE_EXTENSION:
                    file_worker(file)
                    break
        elif _async:
            with Pool(cpu_count(), initializer=async_init) as p:
                p.map(file_worker, self.files)
        else:
            for file in self.files:
                file_worker(file)

        print(f"import done in {datetime.now() - timer}")

def file_worker(file):
    try:
        with transaction.atomic():
            if file.is_file() and file.suffix == SAMPLE_EXTENSION:
                base_name = file.name
                if not EnteredData.objects.filter(
                    type="sample", name=base_name
                ).exists():
                    sample_import = None
                    try:
                        sample_import = SampleImport(file)
                        sample_import.write_to_db()
                        sample_import.move_files(success=True)
                    except Exception as e:
                        if sample_import:
                            sample_import.move_files(success=False)
                        print(f"Error importing file: {base_name}: \n\t {e}\n")
                        raise e
                    else:
                        EnteredData.objects.create(
                            type="sample", name=base_name, date=datetime.now()
                        )
                else:
                    print("already imported ", file)
    except Exception as e:
        print(f"Error importing file: {file}: \n\t {e}\n")
        traceback.print_exc()
