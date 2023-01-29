import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns

class Plot_data():
    """
    This class contains methods for creating and saving exploratory data visualizations.
    It takes two inputs, a dataframe containing the data to be visualized and an optional background dataframe.

    The following visualizations are created:

    gene_characteristics() - Creates bar plots of chromosome distribution and gene type distribution.
    feature_missingness() - Creates a bar plot of fraction of rows with missing data and a csv file containing the percent of missing values for each feature.
    feature_correlation() - Creates a csv file containing the correlation matrix and a csv file of correlation values for all features.
    background_plots() - If a background dataframe is provided, creates visualizations of chromosome and gene type distribution and KS-test p-values for all features.
    """
    
    def __init__(self, df_eda, background_df_eda):
        self.df_eda = df_eda
        self.background_df_eda = background_df_eda

    def gene_characteristics(self):
        print("Plotting data visualizations for input gene list...")
        gene_chars = self.df_eda[['Gene', 'Chr', 'Gene_type']]
        if (gene_chars['Gene_type'] == 1).all():
            value = gene_chars['Gene_type'].iloc[0]
            print(f"Only one gene type in list: {value}")
            chr_plot = sns.barplot(x=gene_chars.Chr.value_counts().index, y=gene_chars.Chr.value_counts())
            chr_plot.set_xlabel("Chromosome")
            chr_plot.set_ylabel("Gene Count")
            plt.savefig('./exploromics_results/Chromosome_distribution.png', format='png', dpi=300, bbox_inches = "tight")
        else:
            gene_type = gene_chars.groupby(['Gene_type']).count()
            chr_plot = sns.barplot(x=gene_chars.Chr.value_counts().index, y=gene_chars.Chr.value_counts())
            chr_plot.set_xlabel("Chromosome")
            chr_plot.set_ylabel("Gene Count")
            plt.savefig('./exploromics_results/Chromosome_distribution.png', format='png', dpi=300, bbox_inches = "tight")
            gene_type_plot = sns.barplot(x=gene_type.index, y=gene_type.Gene)
            gene_type_plot.set_xlabel("Gene type")
            gene_type_plot.set_ylabel("Gene Count")
            plt.xticks(rotation=45)
            plt.savefig('./exploromics_results/Gene_types.png', format='png', dpi=300, bbox_inches = "tight")
        
    def feature_missingness(self):
        df_eda = self.df_eda
        null_counts = df_eda.isnull().sum() / len(df_eda)
        plt.figure(figsize=(30, 10))
        plt.xticks(np.arange(len(null_counts)) + 0.0, null_counts.index, rotation="vertical")
        plt.ylabel("Fraction of rows with missing data")
        plt.bar(np.arange(len(null_counts)), null_counts, color="dodgerblue")
        plt.savefig('./exploromics_results/Feature_missingness.png', format='png', dpi=300, bbox_inches = "tight")

        percent_missing = df_eda.isnull().sum() * 100 / len(df_eda)
        missing_value_df = pd.DataFrame(
            {"column_name": df_eda.columns, "percent_missing": percent_missing}
        )
        missing_value_df.sort_values("percent_missing", inplace=True)
        missing_value_df.to_csv('./exploromics_results/Feature_missingness.csv', index=False)
      
        
    def feature_correlation(self):
        df_corr = self.df_eda
        df_corr.drop(columns=['Gene_type'], inplace=True, errors='ignore')
        df_corr.set_index("Gene", inplace=True)
        corr = df_corr.corr(method="spearman")
        corr.to_csv('./exploromics_results/correlation_matrix.csv', index=True)
        corr_list = corr.stack().reset_index()
        corr_list.columns = ['Variable_1', 'Variable_2', 'Correlation']
        corr_list.to_csv('./exploromics_results/correlation_list.csv', index=False)


    def background_plots(self):
            if self.background_df_eda is not None:
                if self.background_df_eda.iloc[:, 0].equals(self.df_eda.iloc[:, 0]):
                    pass
                else:
                    print("Plotting data visualizations for background genes...")
                    gene_chars = self.background_df_eda[['Gene', 'Chr', 'Gene_type']]
                    if (gene_chars['Gene_type'] == 1).all():
                        value = gene_chars['Gene_type'].iloc[0]
                        print(f"Only one gene type in list: {value}")
                        chr_plot = sns.barplot(x=gene_chars.Chr.value_counts().index, y=gene_chars.Chr.value_counts())
                        chr_plot.set_xlabel("Chromosome")
                        chr_plot.set_ylabel("Gene Count")
                        plt.savefig('./exploromics_results/background_genes/Background_gene_chromosome_distribution.png', format='png', dpi=300, bbox_inches = "tight")
                    else:
                        gene_type = gene_chars.groupby(['Gene_type']).count()
                        chr_plot = sns.barplot(x=gene_chars.Chr.value_counts().index, y=gene_chars.Chr.value_counts())
                        chr_plot.set_xlabel("Chromosome")
                        chr_plot.set_ylabel("Gene Count")
                        plt.savefig('./exploromics_results/background_genes/Background_gene_chromosome_distribution.png', format='png', dpi=300, bbox_inches = "tight")
                        gene_type_plot = sns.barplot(x=gene_type.index, y=gene_type.Gene)
                        gene_type_plot.set_xlabel("Gene type")
                        gene_type_plot.set_ylabel("Gene Count")
                        plt.xticks(rotation=45)
                        plt.savefig('./exploromics_results/background_genes/Background_gene_types.png', format='png', dpi=300, bbox_inches = "tight")
                    
                    background_df_eda = self.background_df_eda
                    null_counts = background_df_eda.isnull().sum() / len(background_df_eda)
                    plt.figure(figsize=(30, 10))
                    plt.xticks(np.arange(len(null_counts)) + 0.0, null_counts.index, rotation="vertical")
                    plt.ylabel("Fraction of rows with missing data")
                    plt.bar(np.arange(len(null_counts)), null_counts, color="dodgerblue")
                    plt.savefig('./exploromics_results/background_genes/Background_genes_feature_missingness.png', format='png', dpi=300, bbox_inches = "tight")

                    percent_missing = background_df_eda.isnull().sum() * 100 / len(background_df_eda)
                    missing_value_df = pd.DataFrame(
                        {"column_name": background_df_eda.columns, "percent_missing": percent_missing}
                    )
                    missing_value_df.sort_values("percent_missing", inplace=True)
                    missing_value_df.to_csv('./exploromics_results/background_genes/Background_genes_feature_missingness.csv', index=False)
                    df_corr = self.background_df_eda
                    df_corr.drop(columns=['Gene_type'], inplace=True, errors='ignore')
                    df_corr.set_index("Gene", inplace=True)
                    corr = df_corr.corr(method="spearman")
                    corr.to_csv('./exploromics_results/background_genes/Background_genes_correlation_matrix.csv', index=True)
                    corr_list = corr.stack().reset_index()
                    corr_list.columns = ['Variable_1', 'Variable_2', 'Correlation']
                    corr_list.to_csv('./exploromics_results/background_genes/Background_genes_correlation_list.csv', index=False)

        
    def gene_background_distributions(self):
        if self.background_df_eda is not None:
            if self.background_df_eda.iloc[:, 0].equals(self.df_eda.iloc[:, 0]):
                pass
            else:
                self.df_eda.drop(columns=['Chr', 'Gene_type'], inplace=True, errors='ignore')
                self.background_df_eda.drop(columns=['Chr', 'Gene_type'], inplace=True, errors='ignore')
                # Create an empty dataframe to store the p-values
                p_values = pd.DataFrame(columns=['Variable', 'p-value'])
                # Iterate through the columns of the dataframes
                for col in self.df_eda.columns:
                    if col in self.background_df_eda.columns:
                        # Perform the Kolmogorov-Smirnov test
                        stat, p = ks_2samp(self.df_eda[col], self.background_df_eda[col])
                        # Add the p-value to the dataframe
                        p_values = pd.concat([p_values, pd.DataFrame({'Variable': [col], 'p-value': [p]})], ignore_index=True)
                        p_values.to_csv('./exploromics_results/background_genes/Gene_background_Kolmogorov-Smirnov_tests.csv', index=False)


    def run(self):
        self.gene_characteristics()
        self.feature_missingness()
        self.feature_correlation()
        self.background_plots()
        self.gene_background_distributions()

        
if __name__ == '__main__':
    plt_dt = Plot_data()
    plt_dt.run()
                
