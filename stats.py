import argparse
import sys
from pathlib import Path
import uuid
import datetime

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.stats as sp

rv_genes = ["BRCA1", "BRCA2", "CDH1", "PALB2", "TP53"]
neg_control_genes = ["BLM", "CEBPA", "FANCA", "FANCB", "GATA2"]
fraction_results = pd.DataFrame()


def permutation_analysis(gene_list, df_case, df_control, iterations=5000):
    all_genes = df_case.gene
    case_genes_length = len(gene_list)  # e.g. 5 genes

    df_case = df_case.dropna(subset=["gene"]) # drop gene '' and trailing other empty genes (all rows must have gene names)
    df_control = df_control.dropna(subset=["gene"]) # drop gene '' and trailing other empty genes (all rows must have gene names)
    #df_case_std = df_case.std()
    #df_control_std = df_control.std()
    df_fraction = df_case.iloc[:, 1:] / df_control.iloc[:, 1:]
    df_case_mean = df_case.mean()
    df_control_mean = df_control.mean()
    for i in range(1,
                   len(df_case.columns)):  # Go through each frequency column, with 0 being the gene.col and simulate j times.
        j = 0
        column_name = df_case.columns[i]
        vals = []

        while j < iterations:
            indices = np.random.choice(df_fraction.shape[0], case_genes_length, replace=False)
            sampled_rows = df_case[indices]
            total_variants_case = df_case[indices, i]
            total_variants_control = df_control[indices, i]  # get the same rows (variant counts for specific gene) from controls
            if total_variants_case.isna().sum() != case_genes_length and total_variants_control.isna().sum() != case_genes_length:
                # Difference of means
                vals.append((total_variants_control - df_control_mean) / (total_variants_case - df_case_mean))
            else:
                vals.append(np.NaN)
            j += 1
        # sc.monte_carlo_test()
        fraction_results[column_name] = vals
    print(fraction_results)
    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results.columns), 2, sharex=False, tight_layout=True, figsize=(20, 20))

    i = 0
    for frequency in fraction_results.columns:
        q =  df_control[df_control.gene.isin(rv_genes)][frequency].mean() / df_case[df_case.gene.isin(rv_genes)][frequency].mean()
        q_control_group =  df_control[df_control.gene.isin(neg_control_genes)][frequency].mean() / df_case[df_case.gene.isin(neg_control_genes)][frequency].mean()

        fraction_results[frequency + "_log"] = np.log2(fraction_results[frequency].astype(dtype=float).dropna())
        if fraction_results[frequency].any():
            mu, std = sp.norm.fit(fraction_results[frequency].dropna())
            xmin, xmax = fraction_results[frequency + "_log"].dropna().min(), fraction_results[
                frequency + "_log"].dropna().max()
            x = np.linspace(xmin, xmax, 100)
            p = sp.norm.pdf(x, mu, std)

            # axs[i].plot(x, p, 'k', linewidth=2)
            axs[i, 0].hist(fraction_results[frequency + "_log"], density=False,
                           log=False, histtype="stepfilled", stacked=False, bins=200)
            axs[i, 0].set_title(
                '{0}: fit values mu ({1}) std({2}) data_mean={3}'.format(fraction_results[frequency].name + "_log",
                                                                         np.round(mu, 2), np.round(std, 2),
                                                                         np.round(fraction_results[frequency].mean(), 2)))
            axs[i, 0].axvline(np.log2(q), color='orange')
            axs[i, 0].axvline(np.log2(q_control_group), color="red")
            axs[i, 1].hist(fraction_results[frequency], density=False,
                           log=False, histtype="stepfilled", bins=200)
            axs[i, 1].set_title(
                "{0} data_mean {1}".format(fraction_results[frequency].name, np.round(fraction_results[frequency].mean(), 2)))
            axs[i, 1].axvline(q, color='orange')

        i += 1
    plt.show()


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
            out.write()
        sys.stderr.write("Created output file {0}.".format(filename))

    permutation_analysis(pd.DataFrame(rv_genes), rv_df, normal_df)
    end = datetime.datetime.now()
    print("Time: {0}".format(end - start))
