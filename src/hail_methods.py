import json
import os
from pathlib import Path

import hail as hl
import hail.utils
import math
import pandas as pd
from hail.utils import info

from .utils import parse_empty, get_metadata, trim_prefix


def create_related_samples_table(mt: hl.MatrixTable) -> hl.Table:
    hl.utils.info("Reducing input.")
    test_intervals = ['2']
    mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval(x, )
         for x in test_intervals])
    mt = mt.filter_rows(hl.len(mt.alleles) > 2, keep=False)
    hl.utils.info("Starting relatedness analysis using KING!")
    mt = mt.unfilter_entries()
    king_table = hl.king(mt.GT)
    #duplicates = king_table.filter_entries(king_table.phi > 0.4, keep=True)  # Find samples that are similar
    #duplicates.key_cols_by()
    #duplicates.show()
    hl.utils.info("Completed KING!")
    return king_table

def multi_way_union_mts(mts: list, tmp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """Joins MatrixTables in the provided list
    :param list mts: list of MatrixTables to join together
    :param str tmp_dir: path to temporary directory for intermediate results
    :param int chunk_size: number of MatrixTables to join per chunk
    :return: joined MatrixTable
    :rtype: MatrixTable
    """
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []
        for i in range(n_jobs):
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols", "__rows")
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols)))
            merged = merged.annotate_rows(
                __rows=hl.flatten(merged.__rows.map(lambda x: x.__rows)))

            print(merged.aggregate((hl.agg.stats(hl.len(merged.__entries)), hl.len(merged.__cols))))
            next_stage.append(
                merged.checkpoint(
                    os.path.join(tmp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"done stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )

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

def merge_matrix_tables_rows(matrix_tables, phenotype=None):
    """
    Merge all matrix tables into one big matrix table with the same header and the same set of samples.
    :param matrix_tables: List of MatrixTables.
    :param phenotype If phenotype is set, try match MatrixTables where mt.globals.phenotype matches this field.
    Deprecated/unused.
    :return: Merged MatrixTable.
    """
    pass
    filtered_mts = matrix_tables

    # do some downfiltering to select only important entries for merging, INFO fields will contain the full data anyway
    combined_mt = matrix_tables[0].select_entries(matrix_tables[0].AD,  matrix_tables[0].DP,
                                                  matrix_tables[0].GT, matrix_tables[0].VF, matrix_tables[0].AC)
    combined_mt = combined_mt.select_rows(combined_mt.impact,
                                          combined_mt.gene, combined_mt.HGNC_ID, combined_mt.MAX_AF)
    if len(matrix_tables) > 1: #  If there is only one match, don't combine any other tables.
        for mt in matrix_tables[1:]:
            mt_prefix = mt.cols
            mt = mt.select_entries(mt.AD, mt.DP, mt.GT, mt.VF, mt.AC)
            mt = mt.select_rows(mt.impact, mt.gene, mt.HGNC_ID, mt.MAX_AF)
            combined_mt = combined_mt.union_rows(mt, _check_cols=False)
    return combined_mt


def merge_matrix_tables_cols(matrix_tables):
    """
    Merge all matrix tables into one big matrix table with the same header and the same set of samples.
    TODO: change signature to write after merge if --write is active (default=False)
    :param matrix_tables: List of MatrixTables.
    :return: Merged MatrixTable.
    """

    # do some downfiltering to select only important entries for merging, INFO fields will contain the full data anyway
    combined_mt = matrix_tables[0].select_entries(matrix_tables[0].AD,  matrix_tables[0].DP,
                                                   matrix_tables[0].GT, matrix_tables[0].VF, matrix_tables[0].AC)
    for mt in matrix_tables[1:]:
        mt = mt.select_entries(mt.AD, mt.DP, mt.GT, mt.VF, mt.AC)
        combined_mt = combined_mt.union_cols(mt, row_join_type="outer")
    return combined_mt


def reduce_to_2d_table(mt, phenotype=None):
    """
    Reduce the matrix table to a 2D matrix table with gene and frequency as keys.
    TODO: rename function, returns a 3D table.

    :param phenotype: Phenotype that is filtered. Deprecated.
    :param mt: Input MatrixTable.
    :return: Returns a dict with the structure: impact { gene { Struct(gnomad_1, gnomad_1_5, gnomad_5) } }
    """
    #Group by globals (phenotype), group by genes, aggregate all into hl.gp_dosage() * 2 (number of total alleles)
    #Filter cols (tables) where phenotype matches command input
    #TODO: create an anti-set where mt.phenotype != phenotype and write that out as anti-table
    # (for statistical comparisons)
    if phenotype is not None:
        mt = mt.filter_cols(mt.phenotype.matches(phenotype), keep=True)
    #out = mt.aggregate_cols(hl.struct(modifier=hl.struct()))
    """
    out = mt.group_rows_by(mt.gene).aggregate_entries(
        modifier=hl.struct(
            gnomad_1=hl.agg.filter(
                (mt.MAX_AF < 0.01) & (mt.impact.contains(hl.literal("MODIFIER"))),
                hl.agg.sum(mt.AC)),
            gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05) & (
                mt.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(mt.AC)),
            gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05) & (
                mt.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(mt.AC))),
        low=hl.struct(
            gnomad_1=hl.agg.filter(
                (mt.MAX_AF < 0.01) & (mt.impact.contains(hl.literal("LOW"))),
                hl.agg.sum(mt.AC)),
            gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05) & (
                mt.impact.contains(hl.literal("LOW"))), hl.agg.sum(mt.AC)),
            gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05) & (
                mt.impact.contains(hl.literal("LOW"))), hl.agg.sum(mt.AC))),
        moderate=hl.struct(
            gnomad_1=hl.agg.filter(
                (mt.MAX_AF < 0.01) & (mt.impact.contains(hl.literal("MODERATE"))),
                hl.agg.sum(mt.AC)),
            gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05) & (
                mt.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(mt.AC)),
            gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05) & (
                mt.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(mt.AC))),
        high=hl.struct(
            gnomad_1=hl.agg.filter(
                (mt.MAX_AF < 0.01) & (mt.impact.contains(hl.literal("HIGH"))),
                hl.agg.sum(mt.AC)),
            gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05) & (
                mt.impact.contains(hl.literal("HIGH"))), hl.agg.sum(mt.AC)),
            gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05) & (
                mt.impact.contains(hl.literal("HIGH"))), hl.agg.sum(mt.AC))))

    #return out.result()
    """
    #mt.summarize()
    results_dict = mt.aggregate_entries(hl.agg.group_by(mt.impact, hl.agg.group_by(mt.gene, hl.struct(
        gnomad_1=hl.agg.filter((mt.MAX_AF < 0.01), hl.agg.sum(mt.AC)),
        gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05), hl.agg.sum(mt.AC)),
        gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05), hl.agg.sum(mt.AC))))))
    #print(results_dict)
    return results_dict

def create_frequency_bins(inp, _raw_out=True):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input dict with a nested structs.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (Dataframe).
    """
    # Convert the dictionary to DataFrame
    df = pd.DataFrame.from_dict(inp, orient='index')
    # Transpose the data so that genes are row indexes
    df = df.T

    unzipped_df = pd.concat([df['HIGH'].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix('HIGH.'),
                             df['MODERATE'].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix('MODERATE.'),
                             df['LOW'].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix('LOW.'),
                             df['MODIFIER'].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix('MODIFIER.')],
                            axis=1)
    condition = unzipped_df.columns.str.endswith('.0') # Remove weird excess columns of 0 that come from concat
    unzipped_df = unzipped_df.loc[:, ~condition]

    # Genes that aren't groupable by a certain impact (e.g. no HIGH in any frequency bin, would return NA,
    # Turn these into 0
    unzipped_df = unzipped_df.fillna(0)
    if _raw_out:
        with hl.utils.with_local_temp_file("hail_raw_json_frequencies.txt") as path:
            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            with open(path, "w",encoding='utf-8') as f:
                f.writelines(inp)
            hl.utils.info("Wrote raw output of bins dict to {0}".format(path))
    return unzipped_df
