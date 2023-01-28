import sys
import argparse
import pandas as pd
import os
from pathlib import Path

class Config:
    def __init__(self, genelist_file, features_file=None, disease=None, disease_id=None, background_genes=None):
        self.genelist_file = genelist_file
        genelist = pd.read_table(genelist_file, sep='\t', header=0)
        self.genelist = genelist 
        self.features_file = features_file
        self.disease = disease
        self.disease_id = disease_id
        self.background_genes = background_genes

        try:
            os.mkdir("exploromics_results")
        except FileExistsError:
            pass
        try:
            if background_genes is not None:
                os.mkdir("exploromics_results/background_genes")   
        except FileExistsError:
            pass

        if features_file is not None:
            if isinstance(self.self.features, pd.DataFrame):
                pass
            else:
                self.self.features  = pd.read_table(self.self.features, sep='\t', header=0)
        if disease is not None:
            self.disease = disease.lower()
        if disease_id is not None:
            self.disease_id = ''.join(map(str, disease_id))
        if background_genes is not None:
            if isinstance(self.background_genes, pd.DataFrame):
                pass
            else:
                self.background_genes  = pd.read_table(self.background_genes, sep='\t', header=0)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Multi-omic annotation and analysis of a list of genes')
    parser.add_argument('--genelist', type=str, required=True, help='Path to the genelist file')
    parser.add_argument("--features", type=str, required=False, default=None, help='Path to user features')
    parser.add_argument('--disease', type=str, required=False,  nargs="?", default=None, help='Disease name')
    parser.add_argument("--diseaseID", type=str, required=False,  nargs="?", default=None, help='Disease ID (as used by OpenTargets)')
    parser.add_argument("--background", type=str,  required=False, default=None, help='Path to background genelist file')
    #Parse command-line arguments
    args = parser.parse_args()

    #Assign the values of the arguments to variables
    genelist_file = args.genelist
    features_file = args.features
    disease = args.disease
    disease_id = args.diseaseID
    background_genes = args.background
    cfg = Config(genelist_file, features_file, disease, disease_id, background_genes)
    