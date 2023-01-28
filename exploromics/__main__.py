import argparse
from exploromics.config_class import Config

class Exploromics(Config):

    def run(self):
        from exploromics.modules.annotation import Annotation
        from exploromics.modules.plot_data import Plot_data

        self.annotate = Annotation(genelist_file=self.genelist_file, features_file=self.features_file, disease=self.disease, 
                                   disease_id=self.disease_id, background_genes=self.background_genes)
        self.annotate.run()
        self.plt_dt = Plot_data(self.annotate.df_eda, self.annotate.background_df_eda)
        self.plt_dt.run()


def main():
    parser = argparse.ArgumentParser(description='Multi-omic annotation and analysis of a list of genes')
    parser.add_argument('--genelist', type=str, required=True, help='Path to the genelist file')
    parser.add_argument("--features", type=str, nargs="?", default=None, help='Path to user features')
    parser.add_argument('--disease', type=str, nargs="?", default=None, help='Disease name')
    parser.add_argument("--diseaseID", type=str, nargs="?", default=None, help='Disease ID (as used by OpenTargets)')
    parser.add_argument("--background", type=str, nargs="?", default=None, help='Path to background genelist file')
    #Parse command-line arguments
    args = parser.parse_args()

    #Assign the values of the arguments to variables
    genelist_file = args.genelist
    features_file = args.features
    disease = args.disease
    disease_id = args.diseaseID
    background_genes = args.background
    
    exploromics = Exploromics(genelist_file, features_file, disease, disease_id, background_genes)

    exploromics.run()

if __name__ == "__main__":
    main()
