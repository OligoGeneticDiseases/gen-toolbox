import argparse
import datetime
import json
import sys
import uuid
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sp

def load_gene_config(json_file):
    """Load gene configuration from a JSON file."""
    with open(json_file, "r") as file:
        config = json.load(file)
    return config

def load_genes_from_file(file_path):
    """Load genes from a file."""
    with open(file_path, "r") as file:
        genes = [line.strip() for line in file.readlines()]
    return genes

def load_genes_from_json(file_path):
    """Load genes from a JSON file."""
    with open(file_path, "r") as file:
        data = json.load(file)
    return data["intersect_genes"], data["rv_genes"], data["neg_control_genes"]

def extract_number_from_filename(filename: str) -> int:
    """Extract the number from the given filename."""
    match = re.search(r'\d+', filename)
    return int(match.group()) if match else None

def permutation_analysis(df_case, df_control, intersect_genes, rv_genes, neg_control_genes, combination_length=5, iterations=20000):
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

    # Filter out empty genes and keep only the intersecting genes in both dataframes
    df_case = df_case.dropna(subset=["gene"])
    df_case = df_case[df_case.gene.isin(intersect_genes)]
    df_control = df_control.dropna(subset=["gene"])
    df_control = df_control[df_control.gene.isin(intersect_genes)]

    # df_case_std = df_case.std()
    # df_control_std = df_control.std()
    # df_fraction = df_case.iloc[:, 1:] / case_count / df_control.iloc[:, 1:] / control_count

    # Calculate the mean of the case and control dataframes
    df_case_mean = df_case.mean()
    df_control_mean = df_control.mean()

    fraction_results_2 = pd.DataFrame()
    num_columns = df_case.shape[1]

    #-------------------------------
    # Loop through each frequency column
    for i in range(1, num_columns):
        frequency_column = df_case.columns[i]

        # Vectorized sampling: randomly select indices for all iterations at once
        indices = np.random.choice(df_case.shape[0], (iterations, case_genes_length), replace=False)

        # Calculate the sum of total variants for the case and control groups using the sampled indices
        total_variants_case = df_case.iloc[indices, i].sum(axis=1)
        total_variants_control = df_control.iloc[indices, i].sum(axis=1)

        # Vectorized ratio calculation: calculate the ratio for all iterations at once
        ratio = (total_variants_case / case_count) / (total_variants_control / control_count )

        # Set NaN values for iterations where the sum of total variants is zero for either case or control groups
        ratio[np.logical_or(total_variants_case == 0, total_variants_control == 0)] = np.NaN

        # Store the calculated ratios in the fraction_results_2 dataframe
        fraction_results_2[frequency_column] = ratio

    # Print the fraction results
    print(fraction_results_2)

#-------------plotting----------------------

    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results_2.columns), 2, sharex="none", tight_layout=True, figsize=(20, 24))

    i = 0
    for frequency_column in fraction_results_2.columns:

        # Average sample normalization enrichment ratios for "likely impactful" and "likely non-impactful" genes
        q_avg = np.divide(np.sum(df_case[df_case.gene.isin(rv_genes)][frequency_column]) / case_count,
                          np.sum(df_control[df_control.gene.isin(rv_genes)][frequency_column]) / control_count )
        q_avg_control_group = np.divide(np.sum(df_case[df_case.gene.isin(neg_control_genes)][frequency_column]) / case_count,
                                        np.sum(df_control[df_control.gene.isin(neg_control_genes)][frequency_column]) / control_count )
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
    analyse.add_argument("--gene_config", "-g", type=validate_file, help="Gene configuration JSON file path",
                         required=True)
    analyse.add_argument("--out", "-", type=validate_file, help="Output file path", required=False)
    analyse.add_argument("--iterations", "-n", type=int, help="Total permutation iterations to be ran. ")

    parser.add_argument("--gene-config", dest="gene_config", help="Path to the gene configuration JSON file",
                        required=True)

    args = parser.parse_args()

    if args.command == "Analyse":
        start = datetime.datetime.now()

        config = load_gene_config(args.gene_config)
        intersect_genes = config["intersect_genes"]
        rv_genes = config["rv_genes"]
        neg_control_genes = config["neg_control_genes"]

        normal_df = pd.read_csv(args.input1, sep="\t", header=0)
        rv_df = pd.read_csv(args.input2, sep="\t", header=0)

        outlines = []

        case_count = extract_number_from_filename(args.input1.name)
        control_count = extract_number_from_filename(args.input2.name)

        permutation_analysis(rv_df, normal_df, intersect_genes, rv_genes, neg_control_genes)

        if args.out is None and len(outlines) > 0:
            filename = str(uuid.uuid4())
            with open(filename, "w+") as out:
                out.write(outlines)
            sys.stderr.write("Created output file {0}.".format(filename))

        end = datetime.datetime.now()
        print("Time: {0}".format(end - start))

        exit()