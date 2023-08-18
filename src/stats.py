import argparse
import datetime
import json
import os.path
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

def permutation_analysis(df_case, df_control, rv_genes, case_count, control_count, neg_control_genes, intersect_genes=None, combination_length=5, iterations=20000):
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

    # Take only the intersect of two dataframes based on gene column
    df_case = df_case[df_case['gene'].isin(df_control['gene'])].copy()
    df_control = df_control[df_control['gene'].isin(df_control['gene'])].copy()  # Normalize 0 values.


    # This section will try to normalise values so that very small values (1..10) when compared to very large values
    # (100000) would not show fractions < 1. In addition, we want to avoid division by 0.
    # Finally normalise for the cohort size i.e. 1000 samples vs 1500 samples etc.
    # TODO: Figure out if there are any extreme ratios that result from such normalisation

    df_case.reset_index(drop=True, inplace=True)
    df_control.reset_index(drop=True, inplace=True)
    columns_to_add = df_case.columns[1:]
    df_case[columns_to_add] = df_case[columns_to_add].multiply(1000)
    df_case[columns_to_add] = df_case[columns_to_add].add(0.00001)
    df_case[columns_to_add] = df_case[columns_to_add].divide(case_count)
    df_control[columns_to_add] = df_control[columns_to_add].multiply(1000)
    df_control[columns_to_add] = df_control[columns_to_add].add(0.00001)
    df_control[columns_to_add] = df_control[columns_to_add].divide(control_count)

    if intersect_genes is not None:
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

    fraction_results_2 = pd.DataFrame(columns=df_case.columns[1:].tolist())
    num_columns = df_case.shape[1]

    #-------------------------------
    # Loop through each frequency column
    # TODO: This secion has been reverted from the vectorized solution to a loop, try and vectorize
    for j in range(1, num_columns):
        frequency_column = df_case.columns[j]
        ratios = list()
        for i in range(1, iterations):

            # Vectorized sampling: randomly select indices for all iterations at once
            indices = np.random.choice(df_case.index, size=case_genes_length, replace=False)
            # Calculate the sum of total variants for the case and control groups using the sampled indices
            total_variants_case = df_case.iloc[indices, j].sum()
            total_variants_control = df_control.iloc[indices, j].sum()

            ratio = np.NaN
            if total_variants_control >= 1 and total_variants_case >= 1:
                ratio = total_variants_case / total_variants_control
            ratios.append(ratio)
            # Set NaN values for iterations where the sum of total variants is zero for either case or control groups
            #ratio[np.logical_or(total_variants_case == 0, total_variants_control == 0)] = np.NaN

            # Store the calculated ratios in the fraction_results_2 dataframe
        fraction_results_2[frequency_column] = ratios

        # Print the fraction results
    print(fraction_results_2)
    print("Expected ratio cases / controls: {0}".format(case_count / control_count))

#-------------plotting----------------------
# This section will plot the results, removed the dual-plot graph, it did not show any new relevant information.

    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results_2.columns), 1, sharex="none", tight_layout=True, figsize=(12, 24))

    i = 0
    for frequency_column in fraction_results_2.columns:

        # Average sample normalization enrichment ratios for "likely impactful" and "likely non-impactful" genes
        q_avg = np.divide(np.sum(df_case[df_case.gene.isin(rv_genes)][frequency_column]),
                          np.sum(df_control[df_control.gene.isin(rv_genes)][frequency_column]) )
        q_avg_control_group = np.divide(np.sum(df_case[df_case.gene.isin(neg_control_genes)][frequency_column]),
                                        np.sum(df_control[df_control.gene.isin(neg_control_genes)][frequency_column]))
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
            axs[i].hist(fraction_results_2[frequency_column + "_log2"], density=False,
                           log=False, histtype="stepfilled", stacked=False, bins=500)
            axs[i].set_title(
                '{0}: data_mean={1} case={2} control={3}'.format(
                    fraction_results_2[frequency_column].name + "_log2",
                    np.round(fraction_results_2[frequency_column + "_log2"].mean(), 4),
                    np.round(np.log2(q_avg), 2), np.round(np.log2(q_avg_control_group), 2)))
            axs[i].axvline(np.log2(q_avg), color='orange')
            axs[i].axvline(np.log2(q_avg_control_group), color="red")

        i += 1
    fig.text(s="Positive set: {0}".format(rv_genes), x=0.2, y=0.5, color="orange", ha="left")
    fig.text(s="Negative set: {0}".format(neg_control_genes), x=0.2, y=0.005, color="red", ha="left")

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

    args = parser.parse_args()

    if args.command == "Analyse":
        start = datetime.datetime.now()

        config = load_gene_config(args.gene_config)
        intersect_genes = config["intersect_genes"]
        rv_genes = config["rv_genes"]
        neg_control_genes = config["neg_control_genes"]

        normal_df = pd.read_csv(args.input1, sep=",", header=0)
        normal_df.columns.values[0] = "gene"
        rv_df = pd.read_csv(args.input2, sep=",", header=0)
        rv_df.columns.values[0] = "gene"
        outlines = []

        case_count = extract_number_from_filename(args.input1.name)
        control_count = extract_number_from_filename(args.input2.name)

        permutation_analysis(rv_df, normal_df, case_count=case_count, control_count=control_count, iterations=args.iterations, rv_genes=rv_genes, neg_control_genes=neg_control_genes, combination_length=20)
        # TODO: store the output
        if args.out is None and len(outlines) > 0:
            filename = str(uuid.uuid4())
            with open(filename, "w+") as out:
                out.write(outlines)
            sys.stderr.write("Created output file {0}.".format(filename))

        end = datetime.datetime.now()
        print("Time: {0}".format(end - start))

        exit()