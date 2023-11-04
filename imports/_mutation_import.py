import collections
import csv

import os
import sys
import pickle
from typing import Any
from typing import Dict
from typing import Generator
from typing import Iterator
from typing import List
from typing import Optional
from typing import Set
from typing import Union


def read_pickle(fname):
        with open(fname, "rb") as handle:
            return pickle.load(handle, encoding="bytes")
def read_var_file_and_insert_into_Mutation_table(var_file):
    var_row_list = []
    if not var_file is None: #if path to var file in .sample file
        with open(var_file, "r") as handle:
            for line in handle:
                if line == "//":
                    break
                vardat = line.strip("\r\n").split("\t")
                var_row_list.append((
                    vardat[4],  # element id
                    vardat[0],  # ref
                    vardat[3],  # alt
                    vardat[1],  # start
                    vardat[2],  # end
                    vardat[6])  # frameshift
                )
            if line != "//":
                sys.exit(
                    "cache error: corrupted file ("
                    + sample_data["var_file"]
                    + ")"
                )
    #var_row_list for INSERT into Variant table completly based on .var file
    #alignment to variant --> fill ALIGNMENT table based on .sample file "seqhash" and "sourceid", store "alignemnt_id"
    #use alignment_id for every variant in .var file and fill "alignment2variant" table
    return var_row_list
    
#SEQUENCE, ALIGNMENT, MUTATION, ALIGNMENT2MUTATION, SAMPLE tables
def process_sample_dict(sample_dict):
    #1. Fill SEQUENCE table get seqhash_id (check uniqness, else skip)
    seq_row = [sample_dict["seqhash"]]
    #2. Fill SAMPLE table
    sample_row = [seqhash_id, sample_dict['name']]
    #3. get sequence_id and fill ALIGNMENT table
    alignment_row = [seqid, sample_dict["source_id"]]
    # 4 get alignment_id and variant file
    # 5 rows for Mutation table insertion
    var_row_list = read_var_file_and_insert_into_Mutation_table(
        sample_dict['var_file'],
        alignment_id)
    # 6 add columns to rows:
    for row in var_row_list:
        parent_id=""
        row = list(row)
        element_id = int(row[0])
        #get sequence with element_id if type == 'source' from ELEMENT
        selected_ref_seq = seq_dict[element_id]
        start_pos = int(row[3])
        if  start_pos <= 0: #start
            pre_ref = ""
        else:
            pre_ref = selected_ref_seq[start_pos - 1]
        # Insert pre_ref at position 2, patent_id (not filled yet, in the end link from protein mutation to variant mutation)
        updated_list.insert(2, pre_ref)
        updated_list.insert(8, parent_id)
        
    # 7 insert var_row_list into mutation table, get var_ids. #check uniquness, insert if not (bulk insert)
    
    # 8 Fill alignment2variant table with var_ids and alignment_id

    # 9 if not sample_data["seqhash"] is None: paranoid test (reconstruction original fasta possible based on mutation table and ref seq)
     

if __name__ == "__main__":
    sample_dict=read_pickle("/home/jules/Documents/PathogenSonar/MN908947_cache/samples/ZW/ZWY1YTI0ZDliOTI1MmNjOTc3NmI5NDc1ZGMyZWViNzYyZjFjMzBjNg.sample")
    for key, value in sample_dict.items():
        print(f"{key}: {value}")

    sample_dict['var_file'] = "/home/jules/Documents/PathogenSonar/MN908947_cache/var/YW/YWZlMGVkM2M4Yzg4MjM2NmQ4ZDFmZDBiYmMyYzVjZWQ2ZWMyNzhiZTQ5NTlkZWI5ODM1Y2U4ODkwZmQwMmM5ZkA3ZDU2MjFjZDNiM2U0OThkMGMyN2ZjY2E5ZDNkM2M1MTY4YzdmM2QzZjk3NzZmMzAwNWM3MDExYmQ5MDA2OGNh.var"
    var_row_list = read_var_file(sample_dict['var_file'])
    for elem in var_row_list:
        if elem[0]!='1':
            print(elem)
