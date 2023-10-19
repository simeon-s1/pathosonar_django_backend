import pathlib
import time
from datetime import datetime
from multiprocessing import Lock, Manager, Pool, cpu_count
from queue import Queue

import django
from django.db import connections, transaction

from rest_api.data_entry.sample_import import SampleImport
from rest_api.models import EnteredData

SAMPLE_EXTENSION = ".sample"


log_lock = Lock()
debug = False
_async = False


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
                    file_worker((file, None))
                    break
        elif _async:
            with Manager() as m:
                queue = m.Queue()
                mutation_lock = m.Lock()
                args = [(file, queue, mutation_lock) for file in self.files]
                connections.close_all()
                with Pool(cpu_count(), initializer=django.setup) as p:
                    p.map_async(file_worker, args)
                    # while not results.ready():
                    #     self.print_estimate_time_left(queue)
                    #     time.sleep(15)
        else:
            mutation_cache = {}
            for file in self.files:
                file_worker(file, mutation_cache)

        print(f"import done in {datetime.now() - timer}")

    def print_estimate_time_left(self, queue: Queue):
        files_left = self.files_number - queue.qsize()
        if files_done := (self.files_number - files_left) > 0:
            time_elapsed = datetime.now() - self.start_timer
            time_left = (time_elapsed / (self.files_number - files_left)) * files_left
            print(f"Estimated time left: {time_left}.")

    def log_error(self, error):
        print(error, flush=True)
        with log_lock:
            self.error_log_file.write(error)
            self.error_log_file.flush()


def file_worker(file, mutation_cache):
    try:
        with transaction.atomic():
            if file.is_file() and file.suffix == SAMPLE_EXTENSION:
                base_name = file.name
                if not EnteredData.objects.filter(
                    type="sample", name=base_name
                ).exists():
                    sample_import = None
                    try:
                        sample_import = SampleImport(file, mutation_cache)
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
