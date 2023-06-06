import os
from pathlib import Path

import hail as hl
import math
import pandas as pd
from hail.utils import info

from .utils import parse_empty


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
def import_and_annotate_vcf_batch(vcfs, annotate=True):
    batch = []
    for vcf in vcfs:
        batch.append(import_and_annotate_vcf(vcf, annotate))
    return batch


def import_and_annotate_vcf(vcf_path, annotate=True):
    """
    Import a VCF file and annotate it using hail.import_vcf() and hail.VEP().

    :param vcf_path: The path to the VCF file.
    :param annotate: Annotates the input VCF file using VEP (default=True). Set to false to skip annotation
    (if already annotated with VEP)
    :return: Annotated MatrixTable.
    """
    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y"})
    mt = hl.import_vcf(vcf_path.__str__(), reference_genome='GRCh37', contig_recoding=contig_recoding)
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
    # TODO: annotate globals here
    mt = mt.drop(mt.vep) # Drop now duplicated field
    mt = mt.filter_entries(mt.VF >= 0.3, keep=True)  # Remove all not ALT_pos < 0.3 / DP > 20
    mt.filter_entries(mt.DP > 30, keep=True)
    return mt

def merge_matrix_tables(matrix_tables):
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
        #combined_mt = combined_mt.union_cols(mt, row_join_type="outer")
        combined_mt = combined_mt.union_rows(mt, _check_cols=False)
    #hl.utils.info("OLIGO: Merged {0} MatrixTables")
    return combined_mt

def reduce_to_2d_table(mt, phenotype=None):
    """
    Reduce the matrix table to a 2D matrix table with gene and frequency as keys.
    TODO: rename function, returns a 3D table.

    :param phenotype: Phenotype that is filtered.
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
    mt.describe()
    results_dict = mt.aggregate_entries(hl.agg.group_by(mt.impact, hl.agg.group_by(mt.gene, hl.struct(
        gnomad_1=hl.agg.filter((mt.MAX_AF < 0.01), hl.agg.sum(mt.AC)),
        gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05), hl.agg.sum(mt.AC)),
        gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05), hl.agg.sum(mt.AC))))))
    #print(results_dict)
    return results_dict

def create_frequency_bins(inp, num_bins=16):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input dict with a nested structs.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (Dataframe).
    """


    for k_impact, impact in inp.items():
        for k_gene, gene in impact.items():
            for frequency in gene.items():
                print(k_impact, k_gene, frequency)
    table = pd.DataFrame(data=inp)
    #table = [["gnomad_1", "gnomad_1_5", "gnomad_5_100"]] = table["HIGH"].apply(lambda x: pd.Series(x*))
    print(pd.json_normalize(inp))
    return table
