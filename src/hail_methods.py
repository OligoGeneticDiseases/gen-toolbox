import hail as hl
from .utils import parse_empty

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
    else: #get the data from the CSQ string
        mt = mt.annotate_rows(VEP_str=mt.info.CSQ.first().split("\\|"))
        mt = mt.annotate_entries(AC=mt.GT.n_alt_alleles(),
                                 VF=hl.float(mt.AD[1] / mt.DP))

        mt = mt.annotate_rows(impact=mt.VEP_str[0],
                              gene=mt.VEP_str[1],
                              HGNC_ID=hl.int(parse_empty(mt.VEP_str[2])),
                              MAX_AF=hl.float(parse_empty(mt.VEP_str[3])))
    # TODO: annotate globals here
    # TODO: filter variants here according to VF/DP ratios etc
    return mt

def merge_matrix_tables(matrix_tables):
    """
    Merge all matrix tables into one big matrix table with the same header and the same set of samples.

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

    :param phenotype: Phenotype that is filtered.
    :param mt: Input MatrixTable.
    :return: Reduced MatrixTable.
    """
    #Group by globals (phenotype), group by genes, aggregate all into hl.gp_dosage() * 2 (number of total alleles)
    #Filter cols (tables) where phenotype matches command input
    #TODO: create an anti-set where mt.phenotype != phenotype and write that out as anti-table
    # (for statistical comparisons)
    if phenotype is not None:
        mt = mt.filter_cols(mt.phenotype == phenotype)
    out = mt.group_rows_by(mt.gene).aggregate(
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
                mt.impact.contains(hl.literal("HIGH"))), hl.agg.sum(mt.AC)))
    )
    return out

def create_frequency_bins(mt, num_bins=16):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input MatrixTable.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (MatrixTable).
    """
    """
    # Group into 16 bins by IMPACT (4 bins) * MAX_AF (4 bins)
    mt = mt.annotate_rows(freq=hl.float64(mt.n_het*2 + mt.n_hom_var) / (2 * hl.agg.count_where(mt.GT.is_defined())))
    mt = mt.annotate_rows(bin=hl.int(hl.min(num_bins - 1, hl.ceil(mt.freq * num_bins))))
    mt = mt.group_rows_by(mt.gene, mt.bin).aggregate_rows(count=hl.agg.count())
    """

    mt = mt.key_cols_by()
    mt = mt.entries()  # Convert from MatrixTable to Table
    #mt.describe()
    #mt.show()
    return mt
