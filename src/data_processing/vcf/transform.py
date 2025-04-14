import hail as hl
import pandas as pd
import os

# Content: VCF data transformation functions.


def reduce_to_2d_table(mt, phenotype=None):
    """
    Reduce the matrix table to a 2D matrix table with gene and frequency as keys.
    TODO: rename function, returns a 3D table.

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
    print(results_dict)
    return results_dict


def create_frequency_bins(inp, _raw_out=True):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input dict with a nested structs.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (Dataframe).
    """
    if isinstance(inp, dict):
        # Convert the dictionary to DataFrame
        df = pd.DataFrame.from_dict(inp, orient="index")
        # Transpose the data so that genes are row indexes
        df = df.T

        unzipped_df = pd.concat(
            [
                df["HIGH"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("HIGH."),
                df["MODERATE"]
                .apply(pd.Series)
                .astype(pd.UInt64Dtype())
                .add_prefix("MODERATE."),
                df["LOW"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("LOW."),
                df["MODIFIER"]
                .apply(pd.Series)
                .astype(pd.UInt64Dtype())
                .add_prefix("MODIFIER."),
            ],
            axis=1,
        )
        condition = unzipped_df.columns.str.endswith(
            ".0"
        )  # Remove weird excess columns of 0 that come from concat
        unzipped_df = unzipped_df.loc[:, ~condition]

        # Genes that aren't groupable by a certain impact (e.g. no HIGH in any frequency bin, would return NA,
        # Turn these into 0
        unzipped_df = unzipped_df.fillna(0)
        if _raw_out:
            with hl.utils.with_local_temp_file("hail_raw_json_frequencies.txt") as path:
                if not os.path.exists(os.path.dirname(path)):
                    os.makedirs(os.path.dirname(path))
                with open(path, "w", encoding="utf-8") as f:
                    f.writelines(inp)
                hl.utils.info("Wrote raw output of bins dict to {0}".format(path))
        return unzipped_df
    else:
        return inp
