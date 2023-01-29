import unittest
import pandas as pd
from exploromics.config_class import Config
from exploromics.annotation import Annotation

class TestAnnotation(unittest.TestCase):
    def setUp(self):
        self.annotation = Annotation(genelist = 'example/genes.txt', disease = 'essential hypertension', disease_id='EFO_0000616', background = 'example/other_genes.txt')
        self.expected_opentargets_df = pd.read_csv('./exploromics/tests/expected_opentargets_df.csv')
        self.expected_annotations_df = pd.read_csv('./exploromics/tests/expected_annotations_df.csv')
        self.expected_background_annotations_df = pd.read_csv('./exploromics/tests/expected_backgronud_annotations_df.csv')
    
    def test_opentargets_association(self):
        opentargets_genes = self.annotation.opentargets_association()
        self.assertIsInstance(opentargets_genes, pd.DataFrame)
        self.assertGreater(opentargets_genes.shape[0], 0)
        self.assertIn('Gene', opentargets_genes.columns)
        self.assertIn('opentargets_association_score', opentargets_genes.columns)
        pd.testing.assert_frame_equal(opentargets_genes, self.expected_opentargets_df)
    
    def test_gene_annotations(self):
        #define the expected columns for the output DataFrame
        expected_columns = ['Gene', 'Chr', 'Gene_length', 'Gene_type', 'GTEx_TPM',
                            'Open_Targets_gene', 'pharmgkb_drugged_gene', 'DGIdb_interaction_group_score',
                            'HIPred_Score', 'SDI_gene', 'pLI', 'O_GLcNAc_Score', 'RVIS_Score',
                            'GWAS_catalog_P-value']

        #call the gene_annotations method and check that it returns a DataFrame
        gene_annotations_df = self.annotation.gene_annotations()
        self.assertIsInstance(gene_annotations_df, pd.DataFrame)

        #check that the output DataFrame has the expected columns
        self.assertEqual(set(gene_annotations_df.columns), set(expected_columns))

        #check that the output DataFrame has at least one row
        self.assertGreater(len(gene_annotations_df), 0)

        #check output matches expected
        pd.testing.assert_frame_equal(gene_annotations_df, self.expected_annotations_df)

    
    def test_background_gene_annotations(self):
        #define the expected columns for the output DataFrame
        expected_columns = ['Gene', 'Chr', 'Gene_length', 'Gene_type', 'GTEx_TPM',
                            'Open_Targets_gene', 'pharmgkb_drugged_gene', 'DGIdb_interaction_group_score',
                            'HIPred_Score', 'SDI_gene', 'pLI', 'O_GLcNAc_Score', 'RVIS_Score',
                            'GWAS_catalog_P-value']

        #call the gene_annotations method and check that it returns a DataFrame
        background_annotations_df = self.annotation.background_annotations()
        self.assertIsInstance(background_annotations_df, pd.DataFrame)

        #check that the output DataFrame has the expected columns
        self.assertEqual(set(background_annotations_df.columns), set(expected_columns))

        #check that the output DataFrame has at least one row
        self.assertGreater(len(background_annotations_df), 0)

        #check output matches expected
        pd.testing.assert_frame_equal(background_annotations_df, self.expected_background_annotations_df)



if __name__ == '__main__':
    unittest.main()