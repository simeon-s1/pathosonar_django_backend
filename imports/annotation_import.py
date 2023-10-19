def get_annotation_ID_by_type(self, effect):
    """
    Get the annotation type id by effect type

    Return id
    """
    try:
        sql = "SELECT id FROM annotation_type WHERE seq_ontology = ?"
        self.cursor.execute(sql, [effect])
        row = self.cursor.fetchone()
        # the mecanism 
        if(row is not None):
            _id = row["id"]
        else:
            # means new combination, we can insert the new data to
            # our database.
            _id = self.insert_effect(effect)
    
    except Exception as e:
        LOGGER.error(e)
        LOGGER.error("Effect keyword: " + effect)
        LOGGER.error(self.cursor.fetchall())
        raise
    return _id

def insert_effect(self,seq_ontology):
    sql = "INSERT IGNORE INTO annotation_type (id, seq_ontology, region) VALUES(?,?,?);"
    self.cursor.execute(sql, [None, seq_ontology, "NONE" ])
    return  self.cursor.lastrowid

def insert_alignment2annotation(self, variant_id, alignment_id, annotation_id):
    sql = "INSERT IGNORE INTO alignment2annotation (variant_id, alignment_id, annotation_id) VALUES(?,?,?);"
    self.cursor.execute(sql, [variant_id, alignment_id, annotation_id])

def read_tsv_snpSift(file_path: str) -> pd.DataFrame:
    """
    Process the TSV file from SnpSift, deduplicate the ANN[*].EFFECT column,
    remove values in ANN[*].IMPACT column, and split the records
    to have one effect per row.
    Returns the modified DataFrame.

    Parameters:
        file_path (str): Path to the input TSV file.

    Returns:
        pd.DataFrame: Modified DataFrame with deduplicated ANN[*].EFFECT column and one effect per row.

    Note:

    """
    try:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file_path, delimiter="\t")
        df = df.drop(["ANN[*].IMPACT"], axis=1, errors="ignore")
        df.rename(columns={"ANN[*].EFFECT": "EFFECT", "ANN[*].ALLELE": "ALT" }, errors="raise", inplace=True)
        # Deduplicate the values in the ANN[*].EFFECT column
        # df["EFFECT"] = df["EFFECT"].str.split(",").apply(set).str.join(",")
        # df['ANN[*].IMPACT'] = ''

        # Split the records into one effect per row
        # df = df.explode('ANN[*].EFFECT')
        df.drop_duplicates(inplace=True)    
        
        # Reset the index
        df = df.reset_index(drop=True)
        # print(df)
        return df
    except KeyError as e:
        LOGGER.error(e)
        LOGGER.error(df.columns)
        raise
    except Exception as e:
        LOGGER.error(e)
        raise
    
def read_sonar_hash(file_path: str):

    with open(file_path, "r") as file:
        data = json.load(file)

    return data


@staticmethod
def process_annotation(db, paired_list, progress=False):
    """

    Steps:
        1. Read annotated txt file and .sonar_hash
        2. Get alignment ID and source element ID
        3. Get variant ID
        4. Insert the 3 IDs into the database.

    Input:
        paired_list = (annotated_file, sonar_hash_file)

    """
    with sonarDBManager(db, readonly=False) as dbm:
        for _tuple in tqdm(
            paired_list,
            desc="Importing annoation data...",
            total=len(paired_list),
            unit="samples",
            bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
            disable=progress,
        ):
            annotated_file, sonar_hash_file = _tuple
            annotated_df = read_tsv_snpSift(annotated_file)
            sonar_hash = read_sonar_hash(sonar_hash_file)
            reference_accession = sonar_hash["reference"]
            sample_dict = sonar_hash["sample_hashes"]
            sample_variant_dict = sonar_hash["sample_variantTable"]

            # Step 2
            for sample_key in sample_dict:
                hash_value = sample_dict[sample_key]
                source_ids_list = dbm.get_element_ids(reference_accession, "source")

                if len(source_ids_list) > 1:
                    LOGGER.error("There is a duplicated element ID!!")
                    sys.exit(1)
                else:
                    source_element_id = source_ids_list[0]

                alnids = dbm.get_alignment_id(hash_value, source_element_id)

                if type(alnids) is list:
                    LOGGER.error(
                        f"Hash value: {hash_value} is not found in the database!!"
                    )
                    sys.exit(1)
                else:
                    alignment_id = alnids

                # Step 3
                sample_variant_list = sample_variant_dict[sample_key]

                # NOTE: MemoryError can be raised if a huge list is converted to a DataFrame
                _df = pd.DataFrame.from_dict(sample_variant_list)

                for row in _df.itertuples():
                    variant_id = getattr(row, "variant_id")

                    selected_var = dbm.get_variant_by_id(variant_id)
                    if selected_var is None:
                        LOGGER.error("No variant was found")
                        LOGGER.warning("This can happen when using a differnet version of database or database instance.")
                        LOGGER.info("Please ensure data import from the corresponding database version.")
                        sys.exit(1)
                    # ref = selected_var["ref"]

                    # VCF: 1-based position
                    # For DEL, we dont do +1
                    ref = (
                        (selected_var["pre_ref"] + selected_var["ref"])
                        if selected_var["alt"] == " "
                        else selected_var["ref"]
                    )

                    if selected_var["alt"] == " ":
                        if selected_var["start"] == 0:
                            start = 1
                        else:
                            start = selected_var["start"]
                    else:
                        start = selected_var["start"] + 1


                    if selected_var["alt"] == " ":
                        if selected_var["start"] == 0:
                            alt = "."
                        else:
                            alt = selected_var["pre_ref"]
                    else:
                        alt = selected_var["alt"]
                    
                    # Handle different kind of SNV (Nucleotide symbol).
                    if alt != "." and len(alt) == 1:
                        alt = sonarDBManager.IUPAC_CODES["nt"][alt.upper()]

                    else:
                        alt = [alt]
                    # Check if it exists in the annotated txt file.
                    selected_rows = annotated_df.loc[
                        (annotated_df["POS"] == start)
                        & (annotated_df["REF"] == ref)
                        & (annotated_df["ALT"].isin(alt))
                    ]
                    # If it does not return any result or more than 1, we should raise an error because
                    # the wrong annotated text file is being used or the database has already been modified.

                    if len(selected_rows) == 0:
                        LOGGER.error(
                            "It appears that the wrong annotated text file is being used "
                            "or the .sonar_hash file is not match to the input "
                            "or the database has already been modified. Please double-check the file "
                            "or database!"
                        )
                        LOGGER.info("Get VAR:")
                        LOGGER.info(selected_var)
                        LOGGER.info("Use for searching a ROW:")
                        LOGGER.info(f"start:{start} , ref:{ref} , alt:{alt}")
                        LOGGER.info("Get DF:")
                        LOGGER.info(f"{annotated_df[annotated_df['POS'] == start]}")
                        sys.exit(1)

                    # Find associated ID from annotationTable.
                    for index, row in selected_rows.iterrows():

                        effect = row["EFFECT"]

                        if effect is None or effect == ".":
                            effect = ""  # Default
                        effect_id = dbm.get_annotation_ID_by_type(effect)

                        # Step 4
                        # Insert into the database
                        dbm.insert_alignment2annotation(
                            variant_id, alignment_id, effect_id
                        )