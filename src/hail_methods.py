import hail as hl

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
    if annotate:
        mt = hl.vep(mt, './vep_settings.json')
        mt = mt.annotate_rows(impact=mt.vep.IMPACT,
                              gene=mt.vep.SYMBOL,
                              HGNC_ID=mt.vep.HGNC_ID,
                              MAX_AF=mt.vep.MAX_AF)
    return mt

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
