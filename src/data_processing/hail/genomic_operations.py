import math
import os

import hail as hl
from hail.utils import info


# Content: Functions specifically related to genomic operations using Hail (from hail_methods.py).


def create_related_samples_table(mt: hl.MatrixTable) -> hl.Table:
    """
    # TODO: Hail/PySpark does not like large sample tables (thousands) and will throw a memory leak
    Try to do KING relatedness inference on a matrixtable containing many samples
    :param mt: A combined matrixtable with entries as samples
    :return: Relatedness table
    """
    hl.utils.info("Reducing input.")
    test_intervals = ["2"]
    mt = hl.filter_intervals(
        mt,
        [
            hl.parse_locus_interval(
                x,
            )
            for x in test_intervals
        ],
    )
    mt = mt.filter_rows(hl.len(mt.alleles) > 2, keep=False)
    hl.utils.info("Starting relatedness analysis using KING!")
    mt = mt.unfilter_entries()
    king_table = hl.king(mt.GT)
    # duplicates = king_table.filter_entries(king_table.phi > 0.4, keep=True)  # Find samples that are similar
    # duplicates.key_cols_by()
    # duplicates.show()
    hl.utils.info("Completed KING!")
    return king_table


def multi_way_union_mts(mts: list, tmp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """
    Joins MatrixTables in the provided list in a multi way fashion, that is tries to preserve the unique structure of
    each MatrixTable. Currently unused.
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
            merged = hl.Table.multi_way_zip_join(
                to_merge, "__entries", "__cols", "__rows"
            )
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
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )
            merged = merged.annotate_rows(
                __rows=hl.flatten(merged.__rows.map(lambda x: x.__rows))
            )

            print(
                merged.aggregate(
                    (hl.agg.stats(hl.len(merged.__entries)), hl.len(merged.__cols))
                )
            )
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
    combined_mt = matrix_tables[0].select_entries(
        matrix_tables[0].AD,
        matrix_tables[0].DP,
        matrix_tables[0].GT,
        matrix_tables[0].VF,
        matrix_tables[0].AC,
    )
    combined_mt = combined_mt.select_rows(
        combined_mt.impact, combined_mt.gene, combined_mt.HGNC_ID, combined_mt.MAX_AF
    )
    if (
        len(matrix_tables) > 1
    ):  #  If there is only one match, don't combine any other tables.
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
    combined_mt = matrix_tables[0].select_entries(
        matrix_tables[0].AD,
        matrix_tables[0].DP,
        matrix_tables[0].GT,
        matrix_tables[0].VF,
        matrix_tables[0].AC,
    )
    for mt in matrix_tables[1:]:
        mt = mt.select_entries(mt.AD, mt.DP, mt.GT, mt.VF, mt.AC)
        combined_mt = combined_mt.union_cols(mt, row_join_type="outer")
    return combined_mt
