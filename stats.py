import argparse
import sys
from pathlib import Path
import uuid
import datetime

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.stats as sp

# rv_genes = ["BRCA1", "BRCA2", "CDH1", "PALB2", "TP53"]
rv_genes = ["BRCA1", "BRCA2", "CHEK2", "PALB2", "ATM"]
# neg_control_genes = ["BLM", "CEBPA", "FANCA", "FANCB", "GATA2"]
neg_control_genes = ["APC", "BMPR1A", "MSH2", "MSH6", "PTEN"]
fraction_results = pd.DataFrame()
fraction_results_2 = pd.DataFrame()

def permutation_analysis(gene_list, df_case, df_control, iterations=50000):
    all_genes = df_case.gene
    case_genes_length = len(gene_list)  # e.g. 5 genes

    df_case = df_case.dropna(
        subset=["gene"])  # drop gene '' and trailing other empty genes (all rows must have gene names)
    df_control = df_control.dropna(
        subset=["gene"])  # drop gene '' and trailing other empty genes (all rows must have gene names)
    # df_case_std = df_case.std()
    # df_control_std = df_control.std()
    df_fraction = df_case.iloc[:, 1:] / 1389 / df_control.iloc[:, 1:] / 826
    df_case_mean = df_case.mean()
    df_control_mean = df_control.mean()
    for i in range(1,
                   len(df_case.columns)):  # Go through each frequency column, with 0 being the gene.col and simulate j times.
        j = 0
        column_name = df_case.columns[i]
        vals = []
        vals2 = []
        while j < iterations:
            # indices = np.random.choice(df_fraction.shape[0], case_genes_length, replace=False)
            # sampled_rows = df_case.iloc[indices]
            # total_variants_case = df_case.iloc[indices, i]
            # total_variants_control = df_control.iloc[indices, i]  # get the same rows (variant counts for specific gene) from controls
            sampled_rows = df_case.sample(case_genes_length, axis=0)
            total_variants_case = sampled_rows.iloc[:case_genes_length, i].sort_index()
            total_variants_control = df_control.iloc[sampled_rows.index - 1, i].sort_index()
            if total_variants_case.sum() != 0 and total_variants_control.sum() != 0:
                # Difference of means
                # Case 5 gene mean / control 5 gene mean (both mean normalized)
                vals.append(np.divide(np.sum(total_variants_case - df_case_mean[i - 1]),
                                      np.sum(total_variants_control - df_control_mean[i - 1])))
                vals2.append(
                   np.divide(np.divide(total_variants_case.sum(), 1389), np.divide(total_variants_control.sum(), 826)))
            else:
                vals.append(np.NaN)
                vals2.append(np.NaN)
            j += 1
        # sc.monte_carlo_test()
        fraction_results[column_name] = vals
        fraction_results_2[column_name] = vals2
    print(fraction_results)
    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results.columns), 2, sharex="none", tight_layout=True, figsize=(20, 24))

    i = 0
    for frequency in fraction_results.columns:
        q = np.sum(df_case[df_case.gene.isin(rv_genes)][frequency] - df_case_mean[frequency]) / np.sum(
            df_control[df_control.gene.isin(rv_genes)][frequency] - df_control_mean[frequency])
        q_control_group = np.sum(
            df_case[df_case.gene.isin(neg_control_genes)][frequency] - df_case_mean[frequency]) / np.sum(
            df_control[df_control.gene.isin(neg_control_genes)][frequency] - df_control_mean[frequency])
        print(frequency, q, q_control_group)
        fraction_results[frequency + "_log"] = np.log2(fraction_results[frequency].dropna())
        if fraction_results[frequency].any():
            mu, std = sp.norm.fit(fraction_results[frequency + "_log"].dropna())
            xmin, xmax = fraction_results[frequency + "_log"].dropna().min(), fraction_results[
                frequency + "_log"].dropna().max()
            x = np.linspace(xmin, xmax, 100)
            p = sp.norm.pdf(x, mu, std)

            # axs[i].plot(x, p, 'k', linewidth=2)
            axs[i, 0].hist(fraction_results[frequency + "_log"], density=True,
                           log=False, histtype="stepfilled", stacked=False, bins=200)
            axs[i, 0].set_title(
                '{0}: fit values mu ({1}) std({2}) data_mean={3} case={4} control={5}'.format(
                    fraction_results[frequency].name + "_log",
                    np.round(mu, 2), np.round(std, 2),
                    np.round(fraction_results[frequency + "_log"].mean(), 2),
                    np.round(q, 2), np.round(q_control_group, 2)))
            axs[i, 0].axvline(np.log2(q), color='orange')
            axs[i, 0].axvline(np.log2(q_control_group), color="red")
            axs[i, 1].hist(fraction_results_2[frequency], density=True,
                           log=True, histtype="stepfilled", bins=200)
            axs[i, 1].set_title(
                "{0} data_mean {1}".format(fraction_results_2[frequency].name,
                                           np.round(fraction_results_2[frequency].mean(), 2)))
            #axs[i, 1].axvline(np.log2(q), color='orange')
            #axs[i, 1].axvline(np.log2(q_control_group), color="red")

        i += 1
    fig.text(s=rv_genes, x=0.2, y=0.05, color="orange", ha="left")
    fig.text(s=neg_control_genes, x=0.7, y=0.05, color="red", ha="left")

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
