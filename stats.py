import argparse
import sys
from pathlib import Path
import uuid
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy.stats import norm


rv_genes = ["BRCA1", "BRCA2", "CDH1", "PALB2", "TP53"]
fraction_results = pd.DataFrame()


def permutation_analysis(gene_list, df_case, df_control, iterations=20000):
    all_genes = df_case.gene
    case_genes_length = len(gene_list)  # e.g. 5 genes
    for i in range(1,
                   len(df_case.columns)):  # Go through each frequency column, with 0 being the gene.col and simulate j times.
        j = 0
        column_name = df_case.columns[i]
        vals = []
        while j < iterations:

            sampled_rows = df_case.sample(case_genes_length, axis=0)
            total_variants_case = sampled_rows.iloc[:case_genes_length, i].sort_index()
            control_samples = df_control.iloc[total_variants_case.index.array, i].sort_index()
            if total_variants_case.isna().sum() != case_genes_length and control_samples.isna().sum() != case_genes_length:
                vals.append((total_variants_case.sum() / control_samples.sum()) / case_genes_length)
            else:
                vals.append(np.NaN)
            j += 1
        # sc.monte_carlo_test()
        fraction_results[column_name] = vals
    print(fraction_results)
    #fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    #plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    #plt.hlines(data=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(len(fraction_results.columns), sharex=True, tight_layout=True, figsize=(20, 45))

    i=0
    for frequency in fraction_results.columns:
        fraction_results[frequency] = np.log2(fraction_results[frequency].astype(dtype=float))
        if fraction_results[frequency].any():
            mu, std = norm.fit(fraction_results[frequency].dropna())
            xmin, xmax = fraction_results[frequency].dropna().min(), fraction_results[frequency].dropna().max()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std)

            #axs[i].plot(x, p, 'k', linewidth=2)
            axs[i].hist(fraction_results[frequency], density=True,
                        log=False, histtype="stepfilled", bins=500)
            axs[i].set_title('{0}: fit values mu ({1}) std({2})'.format(fraction_results[frequency].name, np.round(mu, 2), np.round(std, 2)))
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
