import hail as hl

def import_and_annotate_vcf(vcf_path):
    """
    Import a VCF file and annotate it using hail.import_vcf() and hail.VEP().

    :param vcf_path: The path to the VCF file.
    :return: Annotated MatrixTable.
    """
    mt = hl.import_vcf(vcf_path)
    annotated_mt = hl.vep(mt, './vep_settings.json')
    return annotated_mt

def merge_matrix_tables(matrix_tables):
    """
    Merge all matrix tables into one big matrix table with the same header and the same set of samples.

    :param matrix_tables: List of MatrixTables.
    :return: Merged MatrixTable.
    """
    combined_mt = matrix_tables[0]
    for mt in matrix_tables[1:]:
        combined_mt = combined_mt.union_rows(mt)
    return combined_mt

def reduce_to_2d_table(mt):
    """
    Reduce the matrix table to a 2D matrix table with gene and frequency as keys.

    :param mt: Input MatrixTable.
    :return: Reduced MatrixTable.
    """
    mt = mt.key_rows_by('locus', 'alleles')
    mt = mt.annotate_rows(gene=hl.sorted(mt.vep.transcript_consequences, key=lambda x: x.canonical).gene_symbol)
    mt = mt.group_rows_by(mt.gene).aggregate_rows(n_het=hl.agg.sum(mt.GT.is_het()), n_hom_var=hl.agg.sum(mt.GT.is_hom_var()))
    return mt

def create_frequency_bins(mt, num_bins=16):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input MatrixTable.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (MatrixTable).
    """
    mt = mt.annotate_rows(freq=hl.float64(mt.n_het + 2 * mt.n_hom_var) / (2 * hl.agg.count_where(mt.GT.is_defined())))
    mt = mt.annotate_rows(bin=hl.int(hl.min(num_bins - 1, hl.ceil(mt.freq * num_bins))))
    mt = mt.group_rows_by(mt.gene, mt.bin).aggregate_rows(count=hl.agg.count())
    return mt
