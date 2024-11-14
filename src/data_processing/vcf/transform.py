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
    # TODO: create an anti-set where mt.phenotype != phenotype and write that out as anti-table
    # (for statistical comparisons)
    if phenotype is not None:
        mt = mt.filter_cols(mt.phenotype.matches(phenotype), keep=True)
    # out = mt.aggregate_cols(hl.struct(modifier=hl.struct()))
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
    # mt.summarize()
    if 'impact' in mt.row:
        results_dict = mt.aggregate_entries(
            hl.agg.group_by(
                mt.impact,
                hl.agg.group_by(
                    mt.gene,
                    hl.struct(
                        gnomad_1=hl.agg.filter((mt.MAX_AF < 0.01), hl.agg.sum(mt.AC)),
                        gnomad_1_5=hl.agg.filter((mt.MAX_AF > 0.01) & (mt.MAX_AF < 0.05), hl.agg.sum(mt.AC)),
                        gnomad_5_100=hl.agg.filter((mt.MAX_AF > 0.05), hl.agg.sum(mt.AC)),
                    ),
                ),
            )
        )
    else:
        hl.utils.info("No 'impact' field found. Skipping impact-based grouping.")
        results_dict = {}  # Return an empty dictionary or adjust based on your needs

    return results_dict


def create_frequency_bins(inp, _raw_out=True):
    """
    Create a final output frequency table with the specified number of bins (default: 16).

    :param mt: Input dict with a nested structs.
    :param num_bins: The number of bins to create in the final frequency table (default: 16).
    :return: Final frequency table (Dataframe).
    """
    # Convert the dictionary to DataFrame
    df = pd.DataFrame.from_dict(inp, orient="index")
    # Transpose the data so that genes are row indexes
    df = df.T

    # Attempt to use the original, structured `pd.concat` approach
    try:
        unzipped_df = pd.concat(
            [
                df["HIGH"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("HIGH."),
                df["MODERATE"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("MODERATE."),
                df["LOW"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("LOW."),
                df["MODIFIER"].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix("MODIFIER."),
            ],
            axis=1,
        )
    except KeyError:
        # If a KeyError occurs due to missing impact levels, fall back to a dynamic approach
        data_frames = []
        for impact in ["HIGH", "MODERATE", "LOW", "MODIFIER"]:
            if impact in df:
                impact_df = df[impact].apply(pd.Series).astype(pd.UInt64Dtype()).add_prefix(f"{impact}.")
                data_frames.append(impact_df)
        if data_frames:
            unzipped_df = pd.concat(data_frames, axis=1)
        else:
            unzipped_df = pd.DataFrame()  # Empty DataFrame if no impacts are found

    # Remove any unnecessary columns ending with ".0"
    condition = unzipped_df.columns.str.endswith(".0")
    unzipped_df = unzipped_df.loc[:, ~condition]

    # Fill NaN values with 0 for genes without certain impact levels
    unzipped_df = unzipped_df.fillna(0)

    # Save raw output if specified
    if _raw_out:
        with hl.utils.with_local_temp_file("hail_raw_json_frequencies.txt") as path:
            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            with open(path, "w", encoding="utf-8") as f:
                f.write(str(inp))  # Convert dict to string for writing
            hl.utils.info("Wrote raw output of bins dict to {0}".format(path))

    return unzipped_df
