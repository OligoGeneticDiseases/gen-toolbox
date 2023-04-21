import argparse
import sys
from pathlib import Path
import uuid
import datetime

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.stats as sp

intersect_genes = ["AIP",
                   "APC",
                   "ATM",
                   "BAP1",
                   "BLM",
                   "BMPR1A",
                   "BRCA1",
                   "BRCA2",
                   "BRIP1",
                   "CDC73",
                   "CDH1",
                   "CDK4",
                   "CDKN2A",
                   "CEBPA",
                   "CHEK2",
                   "DDB2",
                   "DICER1",
                   "DIS3L2",
                   "EPCAM",
                   "ERCC2",
                   "ERCC3",
                   "ERCC4",
                   "ERCC5",
                   "FANCA",
                   "FANCB",
                   "FANCC",
                   "FANCD2",
                   "FANCE",
                   "FANCF",
                   "FANCG",
                   "FANCI",
                   "FANCL",
                   "FANCM",
                   "FH",
                   "FLCN",
                   "GATA2",
                   "GPC3",
                   "KIT",
                   "MAX",
                   "MEN1",
                   "MET",
                   "MLH1",
                   "MSH2",
                   "MSH6",
                   "MUTYH",
                   "NBN",
                   "NF1",
                   "NF2",
                   "NSD1",
                   "PALB2",
                   "PHOX2B",
                   "PMS2",
                   "PRKAR1A",
                   "PTCH1",
                   "PTEN",
                   "RAD51C",
                   "RAD51D",
                   "RB1",
                   "RECQL4",
                   "RET",
                   "RHBDF2",
                   "RUNX1",
                   "SDHAF2",
                   "SDHB",
                   "SDHC",
                   "SDHD",
                   "SLX4",
                   "SMAD4",
                   "SMARCB1",
                   "STK11",
                   "SUFU",
                   "TMEM127",
                   "TP53",
                   "TSC1",
                   "TSC2",
                   "VHL",
                   "WT1",
                   "XPA",
                   "XPC",
                   ]
# rv_genes = ["BRCA1", "BRCA2", "CDH1", "PALB2", "TP53"]
rv_genes = ["BRCA1", "BRCA2", "CHEK2", "PALB2", "ATM"]
# neg_control_genes = ["BLM", "CEBPA", "FANCA", "FANCB", "GATA2"]
neg_control_genes = ["APC", "BMPR1A", "MSH2", "MSH6", "PTEN"]

fraction_results_2 = pd.DataFrame()


def permutation_analysis(df_case, df_control, combination_length=5, iterations=20000):
    """
    This permutation analysis is based on Monte Carlo analysis of permutations.
    Sets of N=5 random genes are selected from input variant frequency tables and variant enrichment is compared
    between cases (breast cancer diagnosis) and controls (healthy/family cancer risk).
    Curated gene-sets are placed on the simulation graph. Curated gene sets are selected by current knowledge of gene
    function, ie breast cancer genes are compared against colon cancer genes.
    :param combination_length:
    :param df_case: A numpy table of variants summed in a certain frequency-and-impact bin for each gene, cases.
    :param df_control: A numpy table of variants summed in a certain frequency-and-impact bin for each gene, controls.
    :param iterations: The number of iterations to be simulated across every impact and frequency category (column)
    """
    all_genes = df_case.gene
    case_genes_length = combination_length  # e.g. sets of 5 genes

    df_case = df_case.dropna(
        subset=["gene"])  # drop gene '' and trailing other empty genes (all rows must have gene names)
    df_case = df_case[df_case.gene.isin(intersect_genes)]
    df_control = df_control.dropna(
        subset=["gene"])  # drop gene '' and trailing other empty genes (all rows must have gene names)
    df_control = df_control[df_control.gene.isin(intersect_genes)]
    # df_case_std = df_case.std()
    # df_control_std = df_control.std()
    # df_fraction = df_case.iloc[:, 1:] / 1389 / df_control.iloc[:, 1:] / 826
    df_case_mean = df_case.mean()
    df_control_mean = df_control.mean()
    for i in range(1,
                   len(df_case.columns)):  # Go through each frequency column, with 0 being the gene.col and simulate j times.
        j = 0
        frequency_column = df_case.columns[i]
        vals2 = []
        while j < iterations:
            # indices = np.random.choice(df_fraction.shape[0], case_genes_length, replace=False)
            # sampled_rows = df_case.iloc[indices]
            # total_variants_case = df_case.iloc[indices, i]
            # total_variants_control = df_control.iloc[indices, i]  # get the same rows (variant counts for specific gene) from controls
            sampled_rows = df_case.sample(case_genes_length, axis=0)
            total_variants_case = sampled_rows.iloc[:case_genes_length][frequency_column].sort_index()
            total_variants_control = df_control.loc[sampled_rows.index][frequency_column].sort_index()
            if total_variants_case.sum() != 0 and total_variants_control.sum() != 0:
                # Difference of means
                # Single sample average normalized, averages the variant enrichment level across all samples
                # Depicts the size of the input dataframe (number of cases summed together)
                # Currently a hardcoded value depending on the samples in the case group (1389) vs control group (826)
                vals2.append(
                    np.divide(np.divide(total_variants_case.sum(), 1389), np.divide(total_variants_control.sum(), 826)))
            else:
                vals2.append(np.NaN)
            j += 1
        fraction_results_2[frequency_column] = vals2
    print(fraction_results_2)
    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results_2.columns), 2, sharex="none", tight_layout=True, figsize=(20, 24))

    i = 0
    for frequency_column in fraction_results_2.columns:

        # Average sample normalization enrichment ratios for "likely impactful" and "likely non-impactful" genes
        q_avg = np.divide(np.sum(df_case[df_case.gene.isin(rv_genes)][frequency_column]) / 1389,
                          np.sum(df_control[df_control.gene.isin(rv_genes)][frequency_column]) / 826)
        q_avg_control_group = np.divide(np.sum(df_case[df_case.gene.isin(neg_control_genes)][frequency_column]) / 1389,
                                        np.sum(df_control[df_control.gene.isin(neg_control_genes)][frequency_column]) / 826)
        print("Impact group (av-norm. ): {0}, case_genes_enrichment: {1}, control_genes_enrichment: {2}".format(
            frequency_column,
            q_avg, q_avg_control_group))

        fraction_results_2[frequency_column + "_log2"] = np.log2(fraction_results_2[frequency_column].dropna())
        if fraction_results_2[frequency_column].any():
            #mu, std = sp.norm.fit(fraction_results[frequency + "_log"].dropna())
            xmin, xmax = fraction_results_2[frequency_column + "_log2"].dropna().min(), fraction_results_2[
                frequency_column + "_log2"].dropna().max()
            #x = np.linspace(xmin, xmax, 100)
            #p = sp.norm.pdf(x, mu, std)

            # axs[i].plot(x, p, 'k', linewidth=2)

            # Left sided plot for mean-normalized enrichment ratios
            axs[i, 0].hist(fraction_results_2[frequency_column + "_log2"], density=False,
                           log=False, histtype="stepfilled", stacked=False, bins=500)
            axs[i, 0].set_title(
                '{0}: data_mean={1} case={2} control={3}'.format(
                    fraction_results_2[frequency_column].name + "_log2",
                    np.round(fraction_results_2[frequency_column + "_log2"].mean(), 4),
                    np.round(np.log2(q_avg), 2), np.round(np.log2(q_avg_control_group), 2)))
            axs[i, 0].axvline(np.log2(q_avg), color='orange')
            axs[i, 0].axvline(np.log2(q_avg_control_group), color="red")

            # Right-sided plots for avg-normalized enrichment ratios
            axs[i, 1].hist(fraction_results_2[frequency_column], density=False,
                           log=False, histtype="stepfilled", bins=500)
            axs[i, 1].set_title(
                                '{0}: data_mean={1} case={2} control={3}'.format(
                    fraction_results_2[frequency_column].name + "_log2",
                    np.round(fraction_results_2[frequency_column].mean(), 2),
                    np.round(q_avg, 2), np.round(q_avg_control_group, 2)))
            axs[i, 1].axvline(q_avg, color='orange')
            axs[i, 1].axvline(q_avg_control_group, color="red")

        i += 1
    fig.text(s=rv_genes, x=0.2, y=0.05, color="orange", ha="left")
    fig.text(s=neg_control_genes, x=0.7, y=0.05, color="red", ha="left")

    plt.show()
    print("Done {0} iterations".format(iterations))

def validate_file(arg):
    if (file := Path(arg)).is_file():
        return file
    else:
        raise FileNotFoundError(arg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Statistical analysis command-line tool for gene "
                                          "list and variant burden analysis.")
    subparsers = parser.add_subparsers(title="commands", dest="command")
    analyse = subparsers.add_parser("Analyse", help="Find the Monte Carlo permutation values for a given input.")
    analyse.add_argument("--input1", "-i", type=validate_file, help="Input file path", required=True)
    analyse.add_argument("--input2", "-i2", type=validate_file, help="Input file path", required=True)
    analyse.add_argument("--out", "-", type=validate_file, help="Output file path", required=False)
    analyse.add_argument("--iterations", "-n", type=int, help="Total permutation iterations to be ran. ")
    start = datetime.datetime.now()
    args = parser.parse_args()
    # gene_list = []
    normal_df = pd.read_csv(args.input1, sep="\t", header=0)
    rv_df = pd.read_csv(args.input2, sep="\t", header=0)
    outlines = []
    if args.out is None and len(outlines) > 0:
        filename = str(uuid.uuid4())
        with open(filename, "w+") as out:
            out.write(outlines)
        sys.stderr.write("Created output file {0}.".format(filename))

    permutation_analysis(rv_df, normal_df)
    end = datetime.datetime.now()
    print("Time: {0}".format(end - start))
    exit()