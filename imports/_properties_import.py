import collections
import csv

import os
import sys
from typing import Any
from typing import Dict
from typing import Generator
from typing import Iterator
from typing import List
from typing import Optional
from typing import Set
from typing import Union
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from mpire import WorkerPool
from tqdm import tqdm
import json
from Bio.SeqUtils.CheckSum import seguid
from pathosonar.logging import LoggingConfigurator

def get_database_properties() -> Dict[str, Dict[str, Any]]:
    """
    Returns property data as a dict of dict where key is property name.
    If data is not in the cache, it fetches data from the SQLite database.

    Returns:
        Dict[str, Dict[str, Any]]: A dictionary with property names as keys
        and corresponding property data as values.
    """
    sql = f"SELECT * FROM {self.db_database}.property;"
    self.cursor.execute(sql)
    rows = self.cursor.fetchall()
    database_properties = {} if not rows else {x["name"]: x for x in rows}
    return database_properties

def lineage_sublineage_dict(self) -> Dict[str, str]:
    """
    Property that returns a dictionary mapping lineage to sublineage.
    The dictionary is created based on data read from a SQL query.
    The dictionary is cached for future use, i.e., the SQL query is executed only the first time this property is accessed.

    Returns:
        dict: A dictionary where the keys are lineage and the values are sublineage.
    """
    df = pd.read_sql("SELECT * FROM lineages", self.con)
    lineage_sublineage_dict = dict(zip(df.lineage, df.sublineage))
    return lineage_sublineage_dict

    # FILE HANDLING
@staticmethod
@contextmanager
def open_file_autodetect(file_path: str, mode: str = "r"):
    """
    Opens a file with automatic packing detection.

    Args:
        file_path: The path of the file to open.
        mode: The mode in which to open the file. Default is 'r' (read mode).

    Returns:
        A context manager yielding a file object.
    """
    # Use the magic library to identify the file type
    file_type = magic.from_file(file_path, mime=True)

    if file_type == "application/x-xz":
        file_obj = lzma.open(file_path, mode + "t")  # xz
    elif file_type == "application/gzip":
        file_obj = gzip.open(file_path, mode + "t")  # gz
    elif file_type == "application/zip":
        zip_file = zipfile.ZipFile(file_path, mode)  # zip
        # Assumes there's one file in the ZIP, adjust as necessary
        file_obj = zip_file.open(zip_file.namelist()[0], mode)
    elif file_type == "text/plain":  # plain
        file_obj = open(file_path, mode)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    try:
        yield file_obj
    finally:
        file_obj.close()
        if file_type == "application/zip":
            zip_file.close()

@staticmethod
def _link_columns_to_props(
    col_names: List[str], prop_names: Dict[str, str], quiet: bool
) -> Dict[str, str]:
    """
    Link property columns to their corresponding database properties.

    Args:
        fields: List of column names in the metadata file.
        prop_names: Dictionary mapping database property names to column names in the metadata file.
        quiet: Boolean indicating whether to suppress print statements.

    Returns:
        Dictionary linking file columns (values) to database properties (keys).
    """
    links = {}
    props = sorted(prop_names.keys())
    for prop in props:
        prop_name = prop_names[prop]
        c = col_names.count(prop_name)
        if c == 1:
            links[prop] = prop_name
        elif c > 1:
            LOGGER.error(f"'{prop_name}' is not a unique column.")
            sys.exit(1)
    if "sample" not in links:
        LOGGER.error("Missing 'sample' column assignment.")
        sys.exit(1)
    elif len(links) == 1:
        LOGGER.error("The file does not provide any informative column.")
        sys.exit(1)
    if not quiet:
        for prop in props:
            if prop in links:
                LOGGER.info("  " + prop + " <- " + links[prop])
            else:
                LOGGER.info("  " + prop + " missing")
    return links

@staticmethod
def _get_properties_from_db(db: str) -> Set[str]:
    """Get the properties stored in the database."""
    with sonarUtils.connect_to_db(db) as dbm:
        db_properties = set(dbm.properties.keys())
    db_properties.add("sample")
    return db_properties

@staticmethod
def _get_csv_colnames(fname: str, delim: str) -> List[str]:
    """
    Retrieve the column names of a CSV file.

    Args:
        fname: Filename of the CSV file.
        delim: Delimiter used in the CSV file.

    Returns:
        List of column names.
    """
    with open_file_autodetect(fname) as file:
        return file.readline().strip().split(delim)

# property handling
@staticmethod
def get_db_prop_to_column_name_dict(
    db: str,
    prop_links: List[str],
    autolink: bool,
) -> Dict[str, str]:
    """get property names based on user input."""
    db_properties = _get_properties_from_db(db)
    propnames = {x: x for x in db_properties} if autolink else {}

    for link in prop_links:
        if link.count("=") != 1:
            LOGGER.error(
                "'" + link + "' is not a valid column-to-property assignment."
            )
            sys.exit(1)
        prop, col = link.split("=")
        if prop == "SAMPLE":
            prop = "sample"
        if prop not in db_properties:
            LOGGER.error(
                "Sample property '"
                + prop
                + "' is unknown to the selected database. Use list-props to see all valid properties."
            )
            sys.exit(1)
        propnames[prop] = col
    return propnames

# extract properties form csv/tsv files

@staticmethod
def extract_properties_from_csv_and_tsv(
    csv_files: List[str],
    tsv_files: List[str],
    prop_names: Dict[str, str],
    quiet: bool,
) -> Dict:
    """Process the CSV and TSV files."""
    sample_to_properties_dict = collections.defaultdict(dict)
    # check if necessary
    if not csv_files and not tsv_files:
        return sample_to_properties_dict

    # process files
    file_tuples = [(x, ",") for x in csv_files] + [(x, "\t") for x in tsv_files]
    for fname, delim in file_tuples:
        if not quiet:
            LOGGER.info("linking data from" + fname + "...")
        col_names = _get_csv_colnames(fname, delim)
        col_to_prop_links = _link_columns_to_props(
            col_names, prop_names, quiet
        )
        with open_file_autodetect(fname) as handle:
            csvreader = csv.DictReader(handle, delimiter=delim)
            for row in csvreader:
                sample = row[col_to_prop_links["sample"]]
                for x, v in col_to_prop_links.items():
                    if x != "sample":
                        sample_to_properties_dict[sample][x] = row[v]

    return sample_to_properties_dict


def insert_property_into_property_table(
        self,
        name: str,
        datatype: str,
        querytype: str,
        description: str,
        subject: str,
        standard: Optional[str] = None,
        check_name: bool = True,
    ) -> int:
        """
        Adds a new property and returns the property id.

        Args:
            name (str): The name of the property.
            datatype (str): The data type of the property.
            querytype (str): The query type of the property.
            description (str): The description of the property.
            subject (str): The subject of the property.
            standard (Optional[str], optional): The standard of the property. Defaults to None.
            check_name (bool, optional): Whether to check the property name. Defaults to True.

        adds a new property and returns the property id.

        >>> dbm = getfixture('init_writeable_dbm')
        >>> id = dbm.add_property("NEW_PROP", "text", "text", "my new prop stores text information", "sample")

        """
        name = name.upper()
        if name in self.__illegal_properties:
            LOGGER.error(
                "error: '"
                + str(name)
                + "' is reserved and cannot be used as property name"
            )
            sys.exit(1)

        if check_name and not re.match("^[A-Z][A-Z0-9_]+$", name):
            LOGGER.error(
                "Invalid property name (property names have to start with an letter and can contain only letters, numbers and underscores)"
            )
            sys.exit(1)

        if not re.match("^[a-zA-Z0-9_]+$", name):
            sys.exit(
                "error: invalid property name (property names can contain only letters, numbers and underscores)"
            )
        if name in self.properties:
            sys.exit(
                "error: a property named "
                + name
                + " already exists in the given database."
            )
        try:
            sql = f"INSERT INTO {self.db_database}.property (name, datatype, querytype, description, target, standard) VALUES(?, ?, ?, ?, ?, ?);"
            self.cursor.execute(
                sql, [name, datatype, querytype, description, subject, standard]
            )
            self.__properties = False

            pid = self.properties[name]["id"]
            if standard is not None:
                sql = (
                    f"INSERT INTO {self.db_database}.{self.properties[name]['target']}2property (property_id, value_"
                    + self.properties[name]["datatype"]
                    + ", sample_id) SELECT ?, ?, id FROM sample WHERE 1;"
                )
                vals = [pid, standard]
                self.cursor.execute(sql, vals)

        except sqlite3.Error as error:
            LOGGER.error(f"Failed to insert data into sqlite table ({str(error)}).")
            sys.exit(1)
        except Exception as error:
            LOGGER.error(f"Failed to insert data into table ({str(error)}).")
            sys.exit(1)

        return pid


def insert_properties_into_sample2property_table(
        properties: Dict[str, Dict[str, str]]
    ) -> None:
    illegal = {
        "GENOMIC_PROFILE",
        "SAMPLE_NAME",
        "PROTEOMIC_PROFILE",
        "FRAMESHIFT_MUTATION",
    }
    for sample_name in properties:
        sample_id = dbm.get_sample_id(sample_name)
        if not sample_id:
            continue
        for property_name, property_value in properties[sample_name].items():
            if property_name in illegal:
                LOGGER.error("This proprty name is reserved and cannot be used.")
                sys.exit(1)

            try:
                sql = (
                    f"INSERT INTO sample2property (sample_id, property_id, value_"
                    + database_properties[property_name]["datatype"]
                    + ") VALUES(?, ?, ?)"
                    + " ON DUPLICATE KEY UPDATE value_"
                    + database_properties[property_name]["datatype"]
                    + "=?"
                )
                # tmp solution for insert empty date
                # Incorrect date value: ''
                if database_properties[property_name]["datatype"] == "date":
                    if property_value == "":
                        property_value = None

                self.cursor.execute(
                    sql,
                    [
                        sample_id,
                        database_properties[property_name]["id"],
                        property_value,
                        property_value,
                    ],
                )

            except Exception as e:
                LOGGER.error(e)
                LOGGER.error(f"[insert property] Sample ID:'{str(sample_id)}' cannot be processed")
                sys.exit("If you need an assistance, please contact us.")


def import_properties():
    #User Input:
    #sonar add-prop --name "LINEAGE_LATEST" --dtype text --qtype text --descr "Lineage info."
    #sonar add-prop --name "SEQUENCING_METHOD" --dtype text --qtype text --descr "SEQ method"
    #sonar add-prop --name "POSTAL_CODE" --dtype zip --qtype zip --descr "zip code"
    #sonar add-prop --name "DATE_OF_SAMPLING" --dtype date --qtype date --descr "sampling date"
    #--cols sample=ID

    #PROPERTY TABLE
    #1. pre-set up properties with standard values 
    reserved_properties = [
        ("IMPORTED", "date", "date", "date sample has been imported to the database", "sample",),
        ("MODIFIED", "date", "date", "date when sample data has been modified lastly","sample",)
    ]
    for prop in reserved_properties:
        dbm.add_property(*prop, "sample", check_name=False)

    predifined_properties = [
        ("SEQUENCING_TECH", "text", "text", "Sequencing technologies"),
        ("PROCESSING_DATE", "date", "date", "Submission/Processing date"),
        ("COUNTRY", "text", "text", "Country where a sample belongs to"),
        ("HOST", "text", "text", "e.g., HUMAN"),
        ("ZIP_CODE", "text", "text", "zip code e.g., 33602"),
        ("LAB", "text", "text", "lab id e.g., 11069"),
        ("LINEAGE", "text", "text", "e.g., BA.2 or B.1.1.7"),
        ("TECHNOLOGY", "text", "text", "e.g., ILLUMINA"),
        ("GENOME_COMPLETENESS", "text", "text", "Genome completeness (e.g., partial or complete)"),
        ("LENGTH", "integer", "numeric", "Genome lenght e.g., 197027"),
        ("COLLECTION_DATE", "date", "date", "Keep a sample collection date"),
    ]
    for prop in predifined_properties:
        dbm.add_property(*prop, "sample")

    #2. add properties given by user: name. dtype. qtype, descr, sample or variant, default_value
    # check for illegal properties, check for invalid letters, check if property exists
    insert_property_into_property_table("LINEAGE_LATEST", "text", "text", "Lineage info.", "sample", None)
    insert_property_into_property_table("SEQUENCING_METHOD", "text", "text", "SEQ method", "sample", None)
    insert_property_into_property_table("POSTAL_CODE", "zip", "zip", "zip code", "sample", None)
    insert_property_into_property_table("DATE_OF_SAMPLING", "date", "date", "Lineage info.", "sample", None)
    
    #SAMPLE2PROPERTY table
    #2. prop_links=args.cols, List[str] List of column to property links (formatted as col=prop) to consider for import.
    # ID = sample name in tsv
    # user input
    prop_links="sample=ID"
    autolink=True
    #3. create property user_prop_name_to_db_prop_name_dict 
    # 1. from db properties {prop:prop}
    # 2. taked user input prop_links, split into property and column_name, check if property is in PROPERTY table keys +"sample
    # 3. return dict prop_names= {prop: col_name}"
    db_prop_to_column_name_dict = get_db_prop_to_column_name_dict(db, prop_links, autolink)

    #4. read all csv/tsv files
    # based on file ending use `,`or '\t` seperator for reading
    # reading column names (header)
    # connect column names to property keys with prop_names dict
    sample_to_properties_dict = extract_properties_from_csv_and_tsv(csv_files, tsv_files, db_prop_to_column_name_dict, quiet=True)
    # sample_to_properties_dict: A dictionary of properties, where the key is a sample name and
    # the value is another dictionary of properties for that sample.

    #5. import only samples with id in database
    # dismiss illegal property names (sys exit)
    database_properties = get_database_properties()

    #6. insert sample_id, property_id, property_value (in correct value coulmn) if key pair exist, update value
    insert_properties_into_sample2property_table(sample_to_properties_dict, database_properties)

    #7. insert MODIFIED AND/OR IMPORTED date property

   