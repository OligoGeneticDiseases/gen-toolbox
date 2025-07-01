import hail as hl
import pandas as pd
import os

# Content: VCF data transformation functions.


def reduce_to_2d_table(mt, phenotype=None):
    """
    Reduce the matrix table to a 2D matrix table with gene and frequency as keys.
    TODO: rename function, returns a dict. Phenotype filtering is not used.

    :param phenotype: Phenotype that is filtered. Deprecated.
    :param mt: Input MatrixTable.
    :return: Returns a dict with the structure: impact { gene { Struct(gnomad_1, gnomad_1_5, gnomad_5) } }
    """
    # Group by globals (phenotype), group by genes, aggregate all into hl.gp_dosage() * 2 (number of total alleles)
    # Filter cols (tables) where phenotype matches command input
    # (for statistical comparisons)
    if phenotype is not None:
        mt = mt.filter_cols(mt.phenotype.matches(phenotype), keep=True)
    # out = mt.aggregate_cols(hl.struct(modifier=hl.struct()))

    out = mt.group_rows_by(mt.gene).aggregate_entries(
        modifier=hl.struct(
            gnomad_001=hl.agg.filter(
                (mt.MAX_AF < 0.0001) & (mt.impact.contains(hl.literal("MODIFIER"))),
                hl.agg.sum(mt.AC)),
            gnomad_001_01=hl.agg.filter((mt.MAX_AF > 0.0001) & (mt.MAX_AF < 0.001) & (
                mt.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(mt.AC)),
            gnomad_01_1=hl.agg.filter((mt.MAX_AF > 0.001) & (
                mt.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(mt.AC))),
        low=hl.struct(
            gnomad_001=hl.agg.filter(
                (mt.MAX_AF < 0.0001) & (mt.impact.contains(hl.literal("LOW"))),
                hl.agg.sum(mt.AC)),
            gnomad_001_01=hl.agg.filter((mt.MAX_AF > 0.0001) & (mt.MAX_AF < 0.001) & (
                mt.impact.contains(hl.literal("LOW"))), hl.agg.sum(mt.AC)),
            gnomad_01_1=hl.agg.filter((mt.MAX_AF > 0.001) & (
                mt.impact.contains(hl.literal("LOW"))), hl.agg.sum(mt.AC))),
        moderate=hl.struct(
            gnomad_001=hl.agg.filter(
                (mt.MAX_AF < 0.0001) & (mt.impact.contains(hl.literal("MODERATE"))),
                hl.agg.sum(mt.AC)),
            gnomad_001_01=hl.agg.filter((mt.MAX_AF > 0.0001) & (mt.MAX_AF < 0.001) & (
                mt.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(mt.AC)),
            gnomad_01_1=hl.agg.filter((mt.MAX_AF > 0.001) & (
                mt.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(mt.AC))),
        high=hl.struct(
            gnomad_001=hl.agg.filter(
                (mt.MAX_AF < 0.0001) & (mt.impact.contains(hl.literal("HIGH"))),
                hl.agg.sum(mt.AC)),
            gnomad_001_01=hl.agg.filter((mt.MAX_AF > 0.0001) & (mt.MAX_AF < 0.001) & (
                mt.impact.contains(hl.literal("HIGH"))), hl.agg.sum(mt.AC)),
            gnomad_01_1=hl.agg.filter((mt.MAX_AF > 0.001) & (
                mt.impact.contains(hl.literal("HIGH"))), hl.agg.sum(mt.AC))))

    #return out.result().entries().to_pandas()

    #mt.describe()
    #mt.summarize()
    results_dict = mt.aggregate_entries(
        hl.agg.group_by(
            mt.impact,
            hl.agg.group_by(
                mt.gene,
                hl.struct(
                    gnomad_001=hl.agg.filter((mt.MAX_AF < 0.0001), hl.agg.sum(mt.AC)),
                    gnomad_001_01=hl.agg.filter(
                        (mt.MAX_AF > 0.0001) & (mt.MAX_AF < 0.001), hl.agg.sum(mt.AC)
                    ),
                    gnomad_01_1=hl.agg.filter((mt.MAX_AF > 0.001) & (mt.MAX_AF < 0.01), hl.agg.sum(mt.AC)),
                ),
            ),
        )
    )
    #print(results_dict)
    return results_dict


def create_frequency_bins(inp, _raw_out=True):
    """
    Create a final output frequency table with the specified number of bins (default: 16). Input is either a
    dict or Dataframe. Function is also used to create raw frequency values per sample.

    :param inp: Input dict with a nested structure.
    :param _raw_out: return the unzipped dataframe
    :return: Final frequency table (Dataframe).
    """
    if isinstance(inp, dict):
        records = {}

        for impact_level, gene_dict in inp.items():
            for gene, struct in gene_dict.items():
                if gene not in records:
                    records[gene] = {}
                for bin_name, value in struct.items():
                    records[gene][f"{impact_level}.{bin_name}"] = value

        # Create DataFrame
        df = pd.DataFrame.from_dict(records, orient='index')
        df = df.fillna(0).astype(int)

        if _raw_out:
            with hl.utils.with_local_temp_file("hail_raw_json_frequencies.txt") as path:
                if not os.path.exists(os.path.dirname(path)):
                    os.makedirs(os.path.dirname(path))
                with open(path, "w", encoding="utf-8") as f:
                    f.writelines(inp)
                hl.utils.info("Wrote raw output of bins dict to {0}".format(path))
        return df
    else:
        return inp
def extract_genotype_dataframe(mt: hl.MatrixTable) -> pd.DataFrame:
    # Make sure you have needed identifiers
    mt = mt.select_entries('AC')  # Or other fields: DP, AD, etc. Derived from GT.n_alt_alleles()
    mt = mt.select_rows('gene')  # optional: if gene info is needed

    # Flatten matrix into a table
    flat_table = mt.entries()

    # Collect to pandas (consider limiting if needed)
    df = flat_table.to_pandas()

    # Create a multi-index or unique ID for variants
    df['variant_id'] = df['locus'].astype(str) + "_" + df['alleles'].astype(str)

    # Pivot into genotype matrix format: rows = variant_id, columns = sample_id
    genotype_df = df.pivot(index='variant_id', columns='s', values='AC')

    return genotype_df