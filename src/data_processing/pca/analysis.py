import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def pca_graphing(pca, pca_locs):
    # PCA analysis code is from https://www.jcchouinard.com/pca-with-python/
    # Input file is converted into .tsv format (whitespaces replaced)

    pca_df = pd.read_table(pca, sep=" ")[0:2]
    locs = pd.read_table(pca_locs)

    pca_df_T = pca_df.set_index('SampleID').T.reset_index()
    pca_df_T.columns = ['genotype_id', 'PC1', 'PC2']
    pca_merged = pca_df_T.merge(locs, how="inner", on='genotype_id')
    sns.set()

    sns.lmplot(
        x='PC1',
        y='PC2',
        data=pca_merged,
        hue='population_code',
        fit_reg=False,
        legend=True
    )

    plt.title('2D PCA Graph', y=0.9)
    plt.figure(figsize=(10, 10))
    plt.show()