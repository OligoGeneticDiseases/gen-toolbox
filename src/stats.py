import argparse
import datetime
import json
import re
import sys
import uuid
from pathlib import Path

import scipy.stats as sp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    match = re.search(r"\d+", filename)
    return int(match.group()) if match else None


def permutation_analysis(
    df_case,
    df_control,
    rv_genes,
    case_count,
    control_count,
    neg_control_genes,
    intersect_genes=None,
    combination_length=5,
    iterations=20000,
):
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

    if intersect_genes is not None:
        # Filter out empty genes and keep only the intersecting genes in both dataframes
        df_case = df_case.dropna(subset=["gene"])
        df_case = df_case[df_case.gene.isin(intersect_genes)]
        df_control = df_control.dropna(subset=["gene"])
        df_control = df_control[df_control.gene.isin(intersect_genes)]
        df_case.reset_index(drop=True, inplace=True)
        df_control.reset_index(drop=True, inplace=True)

    #Check if dataframes are of equal length
    if len(df_case.index) != len(df_control.index):
        print("WARNING: Case dataframe length does not match control dataframe length! The intersect of both dataframes will be analysed.")
    # Take only the intersect of two dataframes based on gene column
    intersection_values = set(df_case['gene']).intersection(df_control['gene'])
    df_case = df_case[df_case["gene"].isin(intersection_values)]
    df_case = df_case.sort_values(by="gene")
    df_control = df_control[df_control["gene"].isin(intersection_values)]
    df_control = df_control.sort_values(by="gene")
    df_case.reset_index(drop=True, inplace=True)
    df_control.reset_index(drop=True, inplace=True)
    # The rows must match in order to do index based math
    assert df_case["gene"].equals(df_control["gene"]), "Case and control dataframe indices do not match!"
    # Normalize 0 values.

    # This section will try to normalise values so that very small values (1..10) when compared to very large values
    # (100000) would not show fractions < 1. In addition, we want to avoid division by 0.
    # Finally normalise for the cohort size i.e. 1000 samples vs 1500 samples etc.
    # TODO: Figure out if there are any extreme ratios that result from such normalisation
    expected_ratio = case_count / control_count
    df_case.reset_index(drop=True, inplace=True)
    df_control.reset_index(drop=True, inplace=True)
    columns_to_add = df_case.columns[1:]
    #df_case[columns_to_add] = df_case[columns_to_add].divide(case_count)
    #df_case[columns_to_add] = df_case[columns_to_add].add(0.001)
    #df_case[columns_to_add] = df_case[columns_to_add].multiply(10000)
    #df_case[columns_to_add] = df_case[columns_to_add].multiply(expected_ratio)
    #df_control[columns_to_add] = df_control[columns_to_add].divide(control_count)
    #df_control[columns_to_add] = df_control[columns_to_add].add(0.001)
    #df_control[columns_to_add] = df_control[columns_to_add].multiply(10000)
    #df_control[columns_to_add] = df_control[columns_to_add].multiply(expected_ratio)


    # df_case_std = df_case.std()
    # df_control_std = df_control.std()
    # df_fraction = df_case.iloc[:, 1:] / case_count / df_control.iloc[:, 1:] / control_count


    # Calculate the mean of the case and control dataframes
    df_case_mean = df_case.mean()
    df_control_mean = df_control.mean()
    fraction_results_1 = pd.DataFrame()
    fraction_results_2 = pd.DataFrame(columns=df_case.columns[1:].tolist())
    num_columns = df_case.shape[1]

    print("Case means \n{0}\nControl means \n{1}".format(df_case_mean, df_control_mean))
    print("Expected ratio cases / controls: {0}, log2 {1}".format(expected_ratio, np.log2(expected_ratio)))
    print("Expected ratio cases / controls by group (log2): \n {0}".format(np.log2(df_case_mean / df_control_mean)))

    # Divisible columns, look at the ratios grossly
    divison_result_gross = df_case[columns_to_add] / df_control[columns_to_add]
    r = divison_result_gross[divison_result_gross[:-1] > expected_ratio].dropna(how="all")
    r["gene"] = df_case["gene"]  ## Do nothing


    # -------------------------------
    # Loop through each frequency column
    # TODO: This secion has been reverted from the vectorized solution to a loop, try and vectorize
    k = 0
    df_list = []
    for j in tqdm(range(1, num_columns)):
        frequency_column = df_case.columns[j]
        ratios = list()
        for i in range(1, iterations):
            # Vectorized sampling: randomly select indices for all iterations at once
            indices = np.random.choice(
                df_case.index, size=case_genes_length, replace=False
            )
            # Calculate the sum of total variants for the case and control groups using the sampled indices
            total_variants_case = df_case.iloc[indices, j].sum()
            total_variants_control = df_control.iloc[indices, j].sum()

            ratio = np.NaN
            #if total_variants_control > combination_length*expected_ratio and total_variants_case > combination_length:
            if total_variants_control > 1\
                    and total_variants_case > 1:
                ratio = np.divide(total_variants_case, total_variants_control)
                if np.log2(ratio/expected_ratio) > 0.7 and total_variants_case > case_genes_length and total_variants_control > case_genes_length: # Arbritary factor of 10 i.e. there are 10 times fewer cases than controls
                    re_df = pd.DataFrame({"burden_event":k, "burden_ratio":ratio,"frequency_bin":frequency_column, "gene":df_case.gene[indices], "case":df_case.iloc[indices, j], "control":df_control.iloc[indices, j], "case_control":df_case.iloc[indices, j]/df_control.iloc[indices, j]})
                    k+=1
                    # Save the interesting ratio combinations for further analysis
                    fraction_results_1 = pd.concat([fraction_results_1, re_df.dropna()], ignore_index=True)
                    #fraction_results_1.reset_index(drop=True)
                ratios.append(ratio/expected_ratio)  # Store the ratio value and the indicies that were impactful for that ratio
            else:
                ratios.append(ratio)

            # Store the calculated ratios in the fraction_results_2 dataframe
        fraction_results_2[frequency_column] = ratios

        # Print the fraction results

    #df = fraction_results_1.groupby(["frequency_bin", "gene"])["case_control", "burden_ratio"].agg(pd.Series.mode, dropna=True).reset_index()
    #df = df[df.case_control > 0]
    #filtered_df.to_pickle("{0}_{1}_{2}.pkl".format(case_count, control_count, iterations))
    fraction_results_1.to_csv("{0}_{1}_{2}.csv".format(case_count, control_count, iterations))
    #pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    #print(filtered_df)


    #plt.scatter(data_df=filtered_df, x=filtered_df["case_control"], y=filtered_df["burden_ratio"])
    #plt.show()
    # -------------PCA----------------------
    # scaler = StandardScaler()
    # df[['case_control', 'burden_ratio']] = scaler.fit_transform(df[['case_control', 'burden_ratio']])
    #
    # # Create a PCA instance
    # pca = PCA(n_components=2)
    #
    # # Fit and transform your data_df
    # principalComponents = pca.fit_transform(df[['case_control', 'burden_ratio']])
    # # Create a DataFrame with the principal components
    # principalDf = pd.DataFrame(data_df=principalComponents, columns=['PC1', 'PC2'])
    #
    # # Concatenate the gene labels to the principal components DataFrame
    # finalDf = pd.concat([df['gene'], principalDf], axis=1)
    #
    # # Plot the PCA results
    # plt.scatter(finalDf['PC1'], finalDf['PC2'])
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')
    # plt.title('PCA Analysis')
    # for i, gene in enumerate(finalDf['gene']):
    #     plt.annotate(gene, (finalDf['PC1'][i], finalDf['PC2'][i]))

    #plt.show()


    # -------------plotting----------------------
    # This section will plot the results, removed the dual-plot graph, it did not show any new relevant information.

    # fraction_results["low.gnomad_5_100"] = np.log2(fraction_results["low.gnomad_5_100"].astype(dtype=float))
    # plt.hist(fraction_results["low.gnomad_5_100"], density=True, log=True, histtype="stepfilled", bins=100)
    # plt.hlines(data_df=fraction_results["low.gnomad_5_100_log"])
    fig, axs = plt.subplots(
        len(fraction_results_2.columns),
        1,
        sharex="none",
        tight_layout=False,
        figsize=(12, 24),
    )
    i = 0
    for frequency_column in fraction_results_2.columns:

        #x = fraction_results_2[frequency_column + "_log2"]
        fraction_results_2[frequency_column].dropna(inplace=True)

        if fraction_results_2[frequency_column].any():
            # Average sample normalization enrichment ratios for "likely impactful" and "likely non-impactful" genes
            q_avg = np.divide(
                np.sum(df_case[df_case.gene.isin(rv_genes)][frequency_column]),
                np.sum(df_control[df_control.gene.isin(rv_genes)][frequency_column]),)
            q_avg_control_group = np.divide(
                np.sum(df_case[df_case.gene.isin(neg_control_genes)][frequency_column]),
                np.sum(
                    df_control[df_control.gene.isin(neg_control_genes)][frequency_column]),)
            print(
                "Impact group (av-norm. ): {0}, case_genes_enrichment: {1}, control_genes_enrichment: {2}".format(
                    frequency_column, q_avg, q_avg_control_group
                )
            )
            fraction_results_2[frequency_column + "_log2"] = np.log2(
                fraction_results_2[frequency_column]).dropna()

            # Unused block
            mu, std = sp.norm.fit(fraction_results_2[frequency_column + "_log2"].dropna())
            percentile_99 = np.percentile(fraction_results_2[frequency_column + "_log2"].dropna(), 99)
            xmin, xmax = (
                fraction_results_2[frequency_column + "_log2"].min(),
                fraction_results_2[frequency_column + "_log2"].dropna().max(),
            )
            x = np.linspace(mu - 3 * std, mu + 3*std, 100)
            # End of block

            # Left sided plot for mean-normalized enrichment ratios
            _, bins, _ = axs[i].hist(
                fraction_results_2[frequency_column+ "_log2"],
                density=False,
                log=False,
                histtype="stepfilled",
                stacked=True,
                bins=500,
            )
            p = sp.norm.pdf(bins, mu, std)
            #axs[i].plot(bins, p, 'k', linewidth=2)
            axs[i].set_title(
                "{0}: mean_enrichment={1} positive_genes={2} negative_genes={3} 99th percentile (purple)={4}".format(
                    fraction_results_2[frequency_column].name + "_log2",
                    np.round(fraction_results_2[frequency_column + "_log2"].mean(), 2),
                    np.round(np.log2(q_avg), 2),
                    np.round(np.log2(q_avg_control_group), 2),
                    np.round(percentile_99, 2)
                    )
            )
        # Commented out positive and negative gene sets
            #axs[i].axvline(np.log2(q_avg), color="orange")
            #axs[i].axvline(np.log2(q_avg_control_group), color="red")
            axs[i].axvline(percentile_99, color="purple")
            i +=1
        # fig.text(
        #     s="Positive set: {0}".format(rv_genes), x=0.2, y=0.5, color="orange", ha="left"
        # )
        # fig.text(
        #     s="Negative set: {0}".format(neg_control_genes),
        #     x=0.2,
        #     y=0.005,
        #     color="red",
        #     ha="left",
        # )

    plt.show()




    # i = 0
    # for frequency_column in fraction_results_2.columns:
    #     filtered_df = df[df["frequency_bin"] == frequency_column]
    #     if filtered_df.shape[0] > 0:
    #         #plt.scatter(data=filtered_df["gene"], x=filtered_df["case_control"], y=filtered_df["burden_ratio"])
    #
    #         ax = filtered_df.plot(x='case_control', y='burden_ratio', kind='scatter', figsize=(5, 5), title=frequency_column, legend="gene")
    #         filtered_df[['case_control', 'burden_ratio', 'gene']].apply(lambda x: ax.text(*x), axis=1)
    #         i+=1
    #         plt.show()

    print("Done {0} iterations".format(iterations))


def validate_file(arg):
    if (file := Path(arg)).is_file():
        return file
    else:
        raise FileNotFoundError(arg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Statistical analysis command-line tool for gene "
        "list and variant burden analysis."
    )
    subparsers = parser.add_subparsers(title="commands", dest="command")
    analyse = subparsers.add_parser(
        "Analyse", help="Find the Monte Carlo permutation values for a given input."
    )
    analyse.add_argument(
        "--case", "-i", type=validate_file, help="Input file path cases", required=True
    )
    analyse.add_argument(
        "--control", "-i2", type=validate_file, help="Input file path control", required=True
    )
    analyse.add_argument(
        "--gene_config",
        "-g",
        type=validate_file,
        help="Gene configuration JSON file path",
        required=True,
    )  # TODO: ask if you would like to proceed with src/config/gene_config.json or different
    analyse.add_argument(
        "--out", "-", type=validate_file, help="Output file path", required=False
    )
    analyse.add_argument(
        "--iterations", "-n", type=int, help="Total permutation iterations to be ran. "
    )

    analyse.add_argument(
        "--csv", action="store_true"
    )

    args = parser.parse_args()

    if args.command == "Analyse":
        start = datetime.datetime.now()

        config = load_gene_config(args.gene_config)
        if args.csv:
            intersect_genes = config["intersect_genes_tso"]
            case_df = pd.read_csv(args.case, sep=",", header=0)
            control_df = pd.read_csv(args.control, sep=",", header=0)
        else:
            intersect_genes = config["intersect_genes_tshc"]
            case_df = pd.read_table(args.case, sep="\t", header=0)
            control_df = pd.read_table(args.control, sep="\t", header=0)
        rv_genes = config["rv_genes"]
        neg_control_genes = config["neg_control_genes"]


        case_df.columns.values[0] = "gene"


        control_df.columns.values[0] = "gene"
        outlines = []

        case_count = extract_number_from_filename(args.case.name)
        control_count = extract_number_from_filename(args.control.name)

        permutation_analysis(
            df_case=case_df,
            df_control=control_df,
            case_count=case_count,
            control_count=control_count,
            iterations=args.iterations,
            rv_genes=rv_genes,
            neg_control_genes=neg_control_genes,
            intersect_genes=intersect_genes,
            combination_length=5,
        )
        # TODO: store the output
        if args.out is None and len(outlines) > 0:
            filename = str(uuid.uuid4())
            with open(filename, "w+") as out:
                out.write(
                    outlines
                )  # TODO: Expected type 'str' (matched generic type 'AnyStr'), got 'list' instead
            sys.stderr.write("Created output file {0}.".format(filename))

        end = datetime.datetime.now()
        print("Time: {0}".format(end - start))

        exit()
