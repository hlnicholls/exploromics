import unittest
import pandas as pd
from exploromics.plot_data import Plot_data
import os

class TestPlotData(unittest.TestCase):
    def setUp(self):
        self.df_eda = pd.DataFrame({'Gene': ['A', 'B', 'C', 'D'],
                                   'Chr': [1, 1, 2, 2],
                                   'Gene_type': [1, 2, 3, 4]})
        self.background_df_eda = pd.DataFrame({'Gene': ['A', 'B', 'C', 'D'],
                                              'Chr': [1, 1, 2, 2],
                                              'Gene_type': [1, 2, 3, 4]})
        self.plot_data = Plot_data(self.df_eda, self.background_df_eda)

    def test_gene_characteristics(self):
        self.plot_data.gene_characteristics()
        self.assertTrue(os.path.exists('./exploromics_results/Chromosome_distribution.png'))
        self.assertTrue(os.path.exists('./exploromics_results/Gene_types.png'))
        os.remove('./exploromics_results/Chromosome_distribution.png')
        os.remove('./exploromics_results/Gene_types.png')

    def test_feature_missingness(self):
        self.plot_data.feature_missingness()
        self.assertTrue(os.path.exists('./exploromics_results/Feature_missingness.png'))
        self.assertTrue(os.path.exists('./exploromics_results/Feature_missingness.csv'))
        os.remove('./exploromics_results/Feature_missingness.png')
        os.remove('./exploromics_results/Feature_missingness.csv')
        
    def test_feature_correlation(self):
        self.plot_data.feature_correlation()
        self.assertTrue(os.path.exists('./exploromics_results/Feature_correlation.csv'))
        os.remove('./exploromics_results/Feature_correlation.csv')
        
    def test_background_plots(self):
        self.plot_data.background_plots()
        self.assertTrue(os.path.exists('./exploromics_results/Background_chromosome_distribution.png'))
        self.assertTrue(os.path.exists('./exploromics_results/Background_gene_types.png'))
        self.assertTrue(os.path.exists('./exploromics_results/Background_ks_test_p_values.csv'))
        os.remove('./exploromics_results/Background_chromosome_distribution.png')
        os.remove('./exploromics_results/Background_gene_types.png')
        os.remove('./exploromics_results/Background_ks_test_p_values.csv')

if __name__ == '__main__':
    unittest.main()