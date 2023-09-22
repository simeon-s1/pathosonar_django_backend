def harmonize_seq(seq: str) -> str:
    """
    Harmonizes the input sequence.

    This function trims leading and trailing white spaces, converts the sequence to upper case and
    replaces all occurrences of "U" with "T". It's usually used to standardize the format of a DNA
    or RNA sequence.

    Args:
        seq (str): The input sequence as a string.

    Returns:
        str: The harmonized sequence.
    """
    try:
        return seq.strip().upper().replace("U", "T")
    except AttributeError as e:
        raise ValueError(
            f"Invalid input, expected a string, got {type(seq).__name__}"
            ) from e

# Reference handling
# add reference and fill REFERNCE, MOLECULE, ELEMENT tables
records = [x for x in iter_genbank(reference_gb)]
_add_reference(dbm, records)

@staticmethod
def _add_reference(dbm: sonarDBManager, records: List[Dict[str, Any]]) -> None:
    """Adds references to the database.

    Args:
        dbm (sonarDbManager): Database manager instance.
        records (List[Dict[str, Any]]): List of reference records.
    """
    # get last ref_id from database
    # add new entry to REFERNCE
    ref_id = dbm.add_reference(
        records[0]["accession"],
        records[0]["description"],
        records[0]["organism"],
        1,
        1,
    )

    for i, record in enumerate(records): # Covid only one record, genomes consisting of more segments with multiple records
        # add molecule insert into MOLECULE table
        mol_id = insert_molecule(dbm, ref_id, i, record)

        # add source
        source_id = add_source_to_element_and_elemparts_table(source_element_dict, molecule_id)

        # add genes
        gene_ids = insert_genes_into_element_table(dbm, mol_id, record["gene"], source_id)
        # add cds
        _add_cds(dbm, mol_id, gene_ids, record["cds"])

@staticmethod
def insert_molecule(
    dbm: sonarDBManager, reference_id: int, i: int, record: Dict[str, Any]
) -> int:
    """Adds reference molecules and elements to the database.

    Args:
        dbm (sonarDbManager): Database manager instance.
        ref_id (int): Reference id.
        i (int): Index of the record.
        record (Dict[str, Any]): Reference record.
    Returns
        int: Molecule ID.
    """
    standard = 1 if i == 0 else 0
    if record["symbol"].strip() == "":
        record["symbol"] = accession
    if standard:
        # we update standard to 1 to put every molecule loaded into _source
        sql = "UPDATE molecule SET standard = ? WHERE reference_id = ? AND standard = 1"
        self.cursor.execute(sql, [0, reference_id])
    sql = "INSERT INTO molecule (id, reference_id, type, accession, symbol, description, segment, length, standard) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);"
    self.cursor.execute(
        sql,
        [
            None,
            reference_id,
            record["moltype"],
            record["accession"],
            record["symbol"],
            record["description"],
            i,
            record["length"],
            standard,
        ],
    )
    sql = "SELECT id FROM molecule WHERE accession = ?"
    self.cursor.execute(sql, [accession])
    mid = self.cursor.fetchone()["id"]
    return mid


@staticmethod
def insert_into_element_and_elemparts_table(element_dict, molecule_id,) -> int:
    """Handles source and inserts it into the database.

    Args:
        dbm (sonarDbManager): Database manager instance.
        mol_id (int): Molecule id.
        source (Dict[str, Any]): Source data.

    Returns:
        int: Source id.
    """
    if symbol.strip() == "":
        symbol = accession
    if standard:
        sql = (
            "UPDATE element SET standard = ? WHERE molecule_id = ? AND standard = 1"
        )
        self.cursor.execute(sql, [0, molecule_id])
    sql = "INSERT INTO element (id, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard, parent_id) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
    if not element_dict["strand"]:
        element_dict["strand"] = 0
    
    self.cursor.execute(
        sql,
        [
            None,
            molecule_id,
            element_dict["type"],
            element_dict["accession"],
            element_dict["symbol"],
            element_dict["description"],
            element_dict["start"],
            element_dict["end"],
            element_dict["strand"],
            element_dict["sequence"],
            element_dict["standard"],
            element_dict["parent_id"]
        ],
    )
    sql = "SELECT id FROM element WHERE accession = ? AND molecule_id =?;"
    self.cursor.execute(sql, [accession, molecule_id])
    element_id = self.cursor.fetchone()["id"]
    if parts is not None:
        for part in parts:
            sql = "INSERT IGNORE INTO elempart (element_id, start, end, strand, base, segment) VALUES(?, ?, ?, ?, ?, ?);"
            self.cursor.execute(sql, [eid] + part)

    if element_dict["sequence"] != dbm.get_sequence(element_id):
        LOGGER.error(
            f"Could not recover sequence of '{source['accession']}' (source) form Genbank file."
        )
        sys.exit(1)

    return element_id

@staticmethod
def insert_genes_into_element_and_elemparts_table(
    mol_id: int,
    gene_elem_list: List[Dict[str, Any]],
) -> Dict[str, int]:
    """Handles genes and inserts them into the database.

    Args:
        dbm (sonarDbManager): Database manager instance.
        mol_id (int): Molecule id.
        genes (List[Dict[str, Any]]): List of genes.
        source_id (int): Source id.

    Returns:
        Dict[str, int]: Dictionary of gene ids with gene accessions as keys.
    """
    gene_id_dict = {}
    for elem in gene_elem_list:
        # QUESTION NOTE: some genomes have the repeated gene name
        # What should we do ???, if we skip these duplicated gene symbol
        # so in Elempart Table will not contain all region.
        if elem["accession"] in gene_ids:
            LOGGER.error(
                f"Mutliple entries for '{elem['accession']}' (gene) in Genbank file."
            )
            # sys.exit(1)
            continue

        insert_into_element_and_elemparts_table(element_dict, molecule_id)

        if elem["sequence"] != dbm.extract_sequence(
            gene_id_dict[elem["accession"]], molecule_id=mol_id
        ):
            
            # print(elem["sequence"])
            LOGGER.debug(f"Gene ID: {gene_id_dict}")
            LOGGER.error(
                f"Could not recover sequence of '{elem['accession']}' (gene) from Genbank file"
            )
            sys.exit(1)

    return gene_id_dict

@staticmethod
def insert_cds_into_element_and_elemparts_table(
    mol_id: int,
    gene_ids: Dict[str, int],
    cds: List[Dict[str, Any]],
    transl_table: Optional[int] = 1,
) -> None:
    """Handles coding sequences (CDS) and inserts them into the database.

    Args:
        dbm (sonarDbManager): Database manager instance.
        mol_id (int): Molecule id.
        gene_ids (Dict[str, int]): Dictionary of gene ids.
        cds (List[Dict[str, Any]]): List of coding sequences.
        transl_table (int, optional): Translation table to use.
    """
    for elem in cds:
        elem['parent_id'] = gene_ids[elem["gene"]]
        cds_id = insert_into_element_and_elemparts_table(element_dict, molecule_id)
        if elem["sequence"] != dbm.extract_sequence(
            cds_id, translation_table=transl_table, molecule_id=mol_id
        ):
            LOGGER.error(
                f"Could not recover sequence of '{elem['accession']}' (cds) from Genbank file"
            )
            sys.exit(1)

# GENBANK PARSING
@staticmethod
def _process_segments(
    feat_location_parts: List[Union[FeatureLocation, CompoundLocation]],
    cds: bool = False,
) -> List[List[int]]:
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
            [int(segment.start), int(segment.end), segment.strand, base, i]
        )
        base += round((segment.end - segment.start - 1) / div, 1)
    return segments


@staticmethod
def _extract_gene_feature(feature: SeqFeature, source_seq: str, source_id) -> dict:
    """
    Extracts the details from the gene feature of a genbank record.

    Args:
        feature (SeqFeature): A Biopython SeqFeature instance.
        source_seq (str): Source sequence.

    Returns:
        dict: A dictionary containing the extracted gene details.

    Raises:
        ValueError: If the feature type is not 'gene' or no qualifier for gene accession or symbol found.
    """
    if feature.type != "gene":
        raise ValueError("The provided feature is not a 'gene' feature.")

    if feature.id != "<unknown id>":
        accession = feature.id
    elif "gene" in feature.qualifiers:
        accession = feature.qualifiers["gene"][0]
    elif "locus_tag" in feature.qualifiers:
        accession = feature.qualifiers["locus_tag"][0]
    else:
        raise ValueError("No qualifier for gene accession found.")

    if "gene" in feature.qualifiers:
        symbol = feature.qualifiers["gene"][0]
    elif "locus_tag" in feature.qualifiers:
        symbol = feature.qualifiers["locus_tag"][0]
    else:
        raise ValueError("No qualifier for gene symbol found.")

    gene_details = {
        "type": "gene",
        "accession": accession,
        "symbol": symbol,
        "description": "",
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "strand": feature.strand,
        "sequence": harmonize_seq(feature.extract(source_seq)),
        "parent_id": source_id,
        "standard":0,
        "parts": _process_segments(feature.location.parts),
    }

    return gene_details

@staticmethod
def _extract_cds_feature(feature) -> dict:
    """
    Extracts the details from the CDS (Coding Sequence) feature of a genbank record.

    Args:
        feature (SeqFeature): A Biopython SeqFeature instance.

    Returns:
        dict: A dictionary containing the extracted CDS details.

    Raises:
        ValueError: If the feature type is not 'CDS'.
    """
    if feature.type != "CDS":
        raise ValueError("The provided feature is not a 'CDS' feature.")

    # for x in ["protein_id", "gene"]:
    #    if x not in feature.qualifiers:
    #         raise ValueError(f"Missing {x} qualifier for cds.")
    if feature.id != "<unknown id>":
        accession = feature.id
    elif "protein_id" in feature.qualifiers:
        accession = feature.qualifiers["protein_id"][0]
    elif "locus_tag" in feature.qualifiers:
        accession = feature.qualifiers["locus_tag"][0]
    else:
        raise ValueError("Missing protein_id/locus_tag qualifier for cds.")


    parts = _process_segments(feature.location.parts, True)
    symbol = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else feature.qualifiers["locus_tag"][0]
    sequence = feature.qualifiers.get("translation", [""])[0]
    description = feature.qualifiers.get("product", [""])[0]

    cds_details = {
        "type": "cds",
        "accession": accession,
        "symbol": symbol,
        "description": description,
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "strand": feature.strand,
        "gene": symbol,
        "sequence": sequence,
        "standard": 0,
        "parent_id": None,       # by gene id later
        "parts": parts,
    }

    if sum([abs(x[1] - x[0]) for x in parts]) % 3 != 0:
        raise ValueError(f"The length of cds '{accession}' is not a multiple of 3.")

    return cds_details

    
@staticmethod
def extract_source_feature(gbk_record: SeqRecord, molecule_dict: Dict, source_element_dict:Dict) -> (Dict, Dict):
    """
    Extract source feature from GenBank record.

    Args:
        gb_record: GenBank record.
        gb_data: Dictionary representing GenBank record.

    Returns:
        The updated GenBank data dictionary.

    Raises:
        ValueError: If SeqRecord does not contain exactly one 'source' feature.
    """
    source_feature = list(filter(lambda x: x.type == "source", gbk_record.features))

    source_feature = [x for x in molecule_dict.features if x.type == "source"]
    if len(source_feature) != 1:
        raise ValueError("Expecting exactly one source feature.")

    source_feature = source_feature[0]
    
    source_sequence = str(harmonize_seq(source_feature.extract(molecule_dict.seq)))
    # one source can consist of multiple segments/genomic regions
    # Each segment is represented as a list of integers [start, end, strand, base, index].
    # where needed?? not in database
    source_parts = _process_segments(source_feature.location.parts)
    source_element_dict.update {
                "type": "source"
                "start": int(source_feature.location.start),
                "end": int(source_feature.location.end),
                "sequence": source_sequence,
                "standard": "1",
                "parts": source_parts, # goes in table Elemparts
            }
    molecule_dict.update(
        {
            "moltype": source_feature.qualifiers.get("mol_type", [""])[0],
            "length": len(source_sequence),
            "segment": source_feature.qualifiers.get("segment", [""])[0],
        }
    )
    return molecule_dict, source_element_dict


def import_gbk_file_into_reference_molecule_element_elempart():
    # parent_id not filled yet, planed as link from protein region to nt region
    parent_id = None

    # 1. iter gene bank file: get info for REFERNCE table, ELEMENT table
    # multiple gb_records possible for segmented genomes
    for i, gbk_record in enumerate(SeqIO.parse(fname, "genbank")): 
        #user input if ref is standard
        standard = 1 if i == 0 else 0

        reference_dict = {}
        molecule_dict =  {}
        element_data = {}
        source_element_dict = {}
        if i==0:
        # REFERNCE
            reference_dict["accession"] = (
                gbk_record.name + "." + str(gbk_record.annotations["sequence_version"])
            )
            reference_dict["description"] = gbk_record.description
            reference_dict["organism"] = gbk_record.annotations["organism"]
        # reference_data["type"]= gbk_record.annotations["mol_type"]
            reference_dict["translation_group_id"] = 1
            reference_dict["standard"] = standard

            ref_id = insert_into_ref_table(reference_dict)
        
        #MOLECULE (same as refernce if only one gbk record, else other segments)
        molecule_dict["accession"] = (
                gbk_record.name + "." + str(gbk_record.annotations["sequence_version"])
            )
        molecule_dict["symbol"] = gb_record.annotations.get("symbol", "") #gb_record???
        molecule_dict["description"] = gbk_record.description
        molecule_dict["segment"] = str(i)
        molecule_dict["length"] = ""
        molecule_dict["standard"] = standard

            #source_element_dict["accession"] = molecule_dict["accession"] 
            #source_element_dict["symbol"] = molecule_dict["accession"] 
        source_element_dict["strand"] = ""
        source_element_dict["description"] = ""
        source_element_dict["type"] = "source"
        source_element_dict["standard"]=1
            #source_element_dict["parent_id"]=None

        # check if only one source feature, else error
        # get and unify full ref seq, add  length, type, segment to ref dict
        # and fill elemnt table info for source (seq, start, end, parts? )

        molecule_dict, source_element_dict = extract_source_feature(gbk_record, molecule_dict, source_element_dict)

        molecule_id = insert_molecule(molecule_dict, ref_id)

        #checks if sequence completly in db, standard=1
        source_id = insert_into_element_and_elemparts_table(source_element_dict, molecule_id)

        #ELEMENT
        element_data["gene"] = []
        element_data["cds"] = []
        

        for feat in gb_record.features:
            if "pseudogene" in feat.qualifiers:
                    continue
            if feat.type == "gene":                
                element_data["gene"].append(
                    # accs, symbol, start, end, strand, seq, desc, parts
                    _extract_gene_feature(
                        feat, source_element_dict["sequence"], source_id
                    )
                )
            elif feat.type == "CDS":
                if "pseudogene" in feat.qualifiers:
                    continue
                element_data["cds"].append(
                    _extract_cds_feature(
                        feat
                        )
                )
        gene_id_dict = insert_genes_into_element_and_elemparts_table(molecule_id, element_data["gene"])
        insert_cds_into_element_and_elemparts_table(
            molecule_id: int,
        gene_id_dict: Dict[str, int],
        element_data["cds"] : List[Dict[str, Any]],
            )
        #insert id, accs, descr, organism, transl_id, standard, updates if only standard changed
        #if multiple gbk_record, only for first one into reference

