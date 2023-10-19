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

"""
Mutation Property table:
We can fill this table based on this code snipped:
fill table mutation properties table based in public enum EffectType listed here:
https://github.com/pcingola/SnpEff/blob/dbf43c09fcaefa34c9af6c21aa46e1a07acbd0d4/src/main/java/org/snpeff/snpEffect/EffectType.java#L319

and region by EffectType getGeneRegion()

or we use imapact and seq_ontology from annotation_features while reading ann.vcf files

"""
"""
1. read ann.vcf
columns: segment_acc | nt pos in segment seq | ref_nt | alt_nt | not needed | not needed | mutation_properties

2. process mutation_properties: 
    2.1. split(',') = one annotation, 
    2.2 split('|') annotation = annotation_features, 
         for every annotation process annotation_features: 
                annotation_features[1] = seq_ontology
                annotation_features[2] = impact (new column in mutation properties table)
                (annotation_features[4] = gene_acc --> not gene acc of mutation (to one annotation different genes can be listed --> genes that are effected by this mutation) if annotation_features[5] == "transcript")
                assign annotation region (region column in properties table) via java script (see above)

3. link to mutation id: use columns [segment_acc | nt pos in segment seq | ref nt | alt nt] for finding gene id
     check in gene table, if "nt pos in segment seq" matches range of gene(s) (gene table start -end)
     if true: get corresponding gene_id(s): for all genes:
            look in mutation table for mutation_id with ref == ref_nt, alt == alt_nt and gene_id==gene_id
            fill alignment2annotation table with mutation_id, alignment_id and annotation_id
           !multiple mutation hits possible but rare
     else: == intergenic variant 
            in old database schema --> id of source element = segment_id
            --> in gene table no full genome entry --> need link to segment id --> add another column "segment" to mutation table? either gene_id OR segment_id must be filled
            search in mutation table for mutation_id with ref == ref_nt, alt == alt_nt and segment_id==segment_id
           !only one mutation hit possible

"""


"""
reading .var files:
current element_id = id to full nucleotide sequence reference (value=1) or to cds_acc protein gene
we can fill the table for nt_mutations correct via:
       *find gene(s) by nt pos (start/end) and add id to muation table, multiple genes --> multiple entries
       *no genes in nt pos: add segment_id of reference (new column segment_id? other ideas?)

for cds mutations we have to wait for Notes new output

"""