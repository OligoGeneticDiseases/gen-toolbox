import hail as hl
import json
import os
from pathlib import Path
import hail.utils
import math
import pandas as pd
from hail.utils import info

from .utils import parse_empty, get_metadata, trim_prefix

#Content: Functions related to reading various file formats.

def load_db_batch(mts):
    batch = []
    for mt_path in mts:
        batch.append(load_mt(mt_path.__str__()))
    return batch

def load_mt(mt_path):
    return hl.read_matrix_table(mt_path)

def import_and_annotate_vcf_batch(vcfs, metadata=None, annotate=True, interval=None):
    batch = []
    for vcf in vcfs:
        batch.append(import_and_annotate_vcf(vcf, metadata=metadata, annotate=annotate, interval=interval))
    return batch

def import_and_annotate_vcf(vcf_path, metadata=None, annotate=True, interval=None):
    """
    Import a VCF file and annotate it using hail.import_vcf() and hail.VEP().

    :param vcf_path: The path to the VCF file.
    :param annotate: Annotates the input VCF file using VEP (default=True). Set to false to skip annotation
    (if already annotated with VEP)
    :param interval hl.LocusInterval to import only certain genomic interval
    :return: Annotated MatrixTable.
    """
    prefix = trim_prefix(vcf_path.stem)
    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y"})
    mt = hl.import_vcf(vcf_path.__str__(), reference_genome='GRCh37', contig_recoding=contig_recoding)
    #mt = hl.filter_alleles(mt, lambda allele, i: hl.is_star(mt.alleles[0], allele))
    #updated_info = mt.info.annotate(AC=mt.new_to_old.map(lambda i: mt.info.AC[i - 1]))
    #ds_result = mt.annotate_rows(info=updated_info)
    #mt = ds_result
    if interval is not None:
        mt = hl.filter_intervals(mt,[hl.parse_locus_interval(x, ) for x in interval])
    mt = mt.filter_rows(mt.alleles[1] != '*') # Filter star alleles as these break VEP
    if annotate:
        mt = hl.vep(mt, './vep_settings.json')
        mt = mt.annotate_rows(impact=mt.vep.IMPACT,
                              gene=mt.vep.SYMBOL,
                              HGNC_ID=mt.vep.HGNC_ID,
                              MAX_AF=mt.vep.MAX_AF)
    else: #get the data from the CSQ string
        mt = mt.annotate_rows(vep=mt.info.CSQ.first().split("\\|")) # Convert CSQ string into the expected VEP output


        mt = mt.annotate_rows(impact=mt.vep[0],
                              gene=mt.vep[1],
                              HGNC_ID=hl.int(parse_empty(mt.vep[2])),
                              MAX_AF=hl.float(parse_empty(mt.vep[3])))
    mt = mt.annotate_entries(AC=mt.GT.n_alt_alleles(),
                             VF=hl.float(mt.AD[1] / mt.DP))
    if metadata is not None:
        phen, mut = metadata.get(prefix, ["NA", "NA"])
        if len(phen) == 0: phen = "NA"
        if len(mut) == 0: mut = "NA"
        mt.annotate_globals(metadata=hl.struct(phenotype=phen, mutation=mut))  # Annotate all rows with corresponding meta
    mt = mt.drop(mt.vep) # Drop now duplicated field
    mt = mt.filter_entries(mt.VF >= 0.3, keep=True)  # Remove all not ALT_pos < 0.3 / DP > 30
    mt.filter_entries(mt.DP > 30, keep=True)
    return mt