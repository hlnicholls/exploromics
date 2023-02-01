from exploromics.config_class import Config
import pandas as pd
import numpy as np
import json
import os
from functools import reduce
from collections import defaultdict
import mygene


class Annotation(Config):
    """
    Annotation class extracts and stores gene annotations from different data sources.
    It annotaties the input genelist and can also annotate an optionally background genelist input

    Attributes:
    disease_id (str): disease id
    open_targets_genes (pd.DataFrame): dataframe of open targets gene associations

    Methods:
    opentargets_association(): Identifies overall OpenTargets disease association scores and converts target ids to gene symbols.
    gene_annotations(): Extracts and stores gene annotations from dgidb, gtex, hipred, sdi, pli, o_glncnac, rvis, and pharmgkb.
    background_annotations(
    """
    
    def opentargets_association(self):
        if self.disease_id is not None:
            print("Identifying OpenTargets disease Overall Association Scores...")
            opentarget_df = pd.DataFrame()
            path_to_json = './exploromics/data/OpenTargets_associationByOverallDirect/'
            for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
                with open(path_to_json + file_name) as json_file:
                    data = json_file.readlines()
                    data = [json.loads(line) for line in data]
                    expectedResult = [d for d in data if d['diseaseId'] in self.disease_id]
                    df_ot_dic = pd.DataFrame(expectedResult)
                    opentarget_df = pd.concat([opentarget_df, df_ot_dic], ignore_index=True)
            mg = mygene.MyGeneInfo()
            mg = mygene.MyGeneInfo()
            gene_convert = mg.querymany(list(opentarget_df['targetId']) , scopes='ensembl.gene', fields='symbol', species='human')
            gene_convert = pd.DataFrame(gene_convert)
            gene_convert.rename(columns={'symbol': 'Gene', '_score': 'opentargets_association_score'}, inplace=True)
            open_targets_genes = gene_convert[['Gene', 'opentargets_association_score']]
            self.open_targets_genes = open_targets_genes
            return self.open_targets_genes 


    def gene_annotations(self):
        dgidb_drugs = pd.read_table('./exploromics/data/dgidb_drug_interactions.tsv', usecols = ['gene_name', 'interaction_group_score'])
        dgidb_drugs.rename(columns={'gene_name': 'Gene', 'interaction_group_score': 'DGIdb_interaction_group_score'}, inplace=True)
        gtex = pd.read_table('./exploromics/data/GTEx_TPM.txt', sep=',')
        gtex.drop(columns = 'Name', inplace= True)
        hipred = pd.read_table('./exploromics/data/HIPred.tsv', sep='\t')
        sdi = pd.read_table('./exploromics/data/SDI_genes.txt', sep=',')
        pli = pd.read_table('./exploromics/data/exac_march16_pLI.txt', sep='\t', usecols = ['gene','pLI'])
        pli.rename(columns={'gene': 'Gene'}, inplace=True)
        o_glncnac = pd.read_table('./exploromics/data/O_GLcNAc_Score.txt', sep=',')
        rvis = pd.read_table('./exploromics/data/genic_intolerance_v3_12Mar16.txt', sep='\t', usecols = ['Gene', 'ALL_0.1%'])
        rvis.rename(columns={'ALL_0.1%': 'RVIS_Score'}, inplace=True)
        pharmgkb_drugs = pd.read_table('./exploromics/data/PharmGKB_relationships.tsv', usecols = ['Entity1_name','Entity1_type','Entity2_name','Association'])
        pharmgkb_drugs.rename(columns={'Entity1_name': 'Gene'}, inplace=True)
        pharmgkb_drugs = pharmgkb_drugs[(pharmgkb_drugs['Entity1_type'] == 'Gene') & (pharmgkb_drugs['Association'] == 'associated')]
        pharmgkb_drugs['Entity2_name'] = pharmgkb_drugs['Entity2_name'].str.lower()
        pharmgkb_disease_drugs = pharmgkb_drugs[(pharmgkb_drugs['Entity2_name'] == self.disease)]
        pharmgkb_disease_drugs = pharmgkb_disease_drugs[['Gene']]
        pharmgkb_disease_drugs['pharmgkb_drugged_gene'] = 1
        gwascatalog = pd.read_table('./exploromics/data/gwas_catalog_v1.0-associations_e108_r2022-12-21.tsv', usecols = ['REPORTED GENE(S)', 'DISEASE/TRAIT', 'P-VALUE'], sep='\t')
        gene_chars = pd.read_table('./exploromics/data/hg19Rel92_AllgeneTypes_0kb.txt', sep='\t', header=None)
        gene_chars.rename(columns={0: 'Chr', 1:'Start', 2:'End', 3:'Gene', 4:'Gene_type', 5:'Source',}, inplace=True)
        gene_chars['Gene_length'] = gene_chars['End'] - gene_chars['Start']
        gene_chars = gene_chars[['Gene', 'Chr', 'Gene_length', 'Gene_type']]
        gene_chars['Chr'] = gene_chars['Chr'].astype('Int64')
        df_gene_chars = pd.merge(self.genelist, gene_chars, on='Gene', how='left')
        gene_length_tmp = df_gene_chars[["Gene", "Chr", "Gene_length", "Gene_type"]]
        gene_length = gene_length_tmp.drop_duplicates(subset=['Gene'])
        
        
        if self.disease is not None and self.disease_id is not None:
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, self.open_targets_genes, pharmgkb_disease_drugs,
                                                                            dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis])
        elif self.disease is not None and self.disease_id is None:
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, pharmgkb_disease_drugs,
                                                                     dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis])
        elif self.disease is None and self.disease_id is not None:
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, self.open_targets_genes, dgidb_drugs,
                                                                    hipred, sdi, pli, o_glncnac, rvis])

        elif self.disease is None and self.disease_id is None:
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, dgidb_drugs,
                                                                    hipred, sdi, pli, o_glncnac, rvis])

        
        if self.features_file is not None and self.disease is not None and self.disease_id is not None:
            self.features_file.iloc[:, 1:] = self.features_file.iloc[:, 1:].applymap(int)
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, self.open_targets_genes, pharmgkb_disease_drugs,
                                                                            dgidb_drugs, hipred, sdi, pli, o_glncnac, self.features_file])
            
        elif self.features_file is not None and self.disease is None and self.disease_id is not None:
            self.features_file.iloc[:, 1:] = self.features_file.iloc[:, 1:].applymap(int)
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, self.open_targets_genes, dgidb_drugs,
                                                                           hipred, sdi, pli, o_glncnac, rvis, self.features_file])

            
        elif self.features_file is not None and self.disease is not None and self.disease_id is None:
            self.features_file.iloc[:, 1:] = self.features_file.iloc[:, 1:].applymap(int)
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, pharmgkb_disease_drugs,
                                                                            dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis, self.features_file])

            
        elif self.features_file is not None and self.disease is None and self.disease_id is None:
            self.features_file.iloc[:, 1:] = self.features_file.iloc[:, 1:].applymap(int)
            df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.genelist, gene_length, gtex, dgidb_drugs,
                                                                           hipred, sdi, pli, o_glncnac, rvis, self.features_file]) 

            
        if self.features_file is not None:
            initial_columns = set(df_eda.columns)
            df_eda_deduplicated = df_eda.drop_duplicates(subset=df_eda.columns, keep=False)
            dropped_columns = initial_columns - set(df_eda_deduplicated.columns)
            if dropped_columns:
                print("The following gene annotation columns were dropped (input features match automated integrated features):", dropped_columns)
            df_eda = df_eda.drop_duplicates(subset=df_eda.columns, keep=False)

        else:
            pass
        
        df_eda = df_eda.drop_duplicates(subset=['Gene'])
        df_eda.to_csv('./exploromics_results/multi-omic_integration_output.csv', index=False) 
        self.df_eda = df_eda 
        gene_length_dict = defaultdict(list)

        nan_rows = df_eda[df_eda.isna().all(axis=1)]
        count = nan_rows.shape[0]
        print(f'{count} input genes have all NaN values across annotations')

        for col in df_eda .columns[4:]:
            corr = df_eda [str(col)].corr(df_eda['Gene_length'], method='spearman')
            gene_length_dict[str(col)].append(corr)
            gene_length_corr = pd.DataFrame.from_dict(gene_length_dict)
            gene_length_corrT = gene_length_corr.T
            gene_length_corrT.rename(columns={0: 'Gene length pearson correlation'}, inplace=True)
        gene_length_corrT.to_csv('./exploromics_results/Gene_length_correlation.csv', index=False)
        print("Multi-omic data integrated for input gene list...")
        return self.df_eda

    def background_annotations(self):
        if self.background_genes is not None:
            if self.background_genes['Gene'].equals(self.genelist['Gene']):
                print('Input gene list and background genes are the same, passing background annotation...')
                pass
            else:
                print('Annotating background genes...')
                dgidb_drugs = pd.read_table('./exploromics/data/dgidb_drug_interactions.tsv', usecols = ['gene_name', 'interaction_group_score'])
                dgidb_drugs.rename(columns={'gene_name': 'Gene', 'interaction_group_score': 'DGIdb_interaction_group_score'}, inplace=True)
                gtex = pd.read_table('./exploromics/data/GTEx_TPM.txt', sep=',')
                gtex.drop(columns = 'Name', inplace= True)
                hipred = pd.read_table('./exploromics/data/HIPred.tsv', sep='\t')
                sdi = pd.read_table('./exploromics/data/SDI_genes.txt', sep=',')
                pli = pd.read_table('./exploromics/data/exac_march16_pLI.txt', sep='\t', usecols = ['gene','pLI'])
                pli.rename(columns={'gene': 'Gene'}, inplace=True)
                o_glncnac = pd.read_table('./exploromics/data/O_GLcNAc_Score.txt', sep=',')
                rvis = pd.read_table('./exploromics/data/genic_intolerance_v3_12Mar16.txt', sep='\t', usecols = ['Gene', 'ALL_0.1%'])
                rvis.rename(columns={'ALL_0.1%': 'RVIS_Score'}, inplace=True)
                pharmgkb_drugs = pd.read_table('./exploromics/data/PharmGKB_relationships.tsv', usecols = ['Entity1_name','Entity1_type','Entity2_name','Association'])
                pharmgkb_drugs.rename(columns={'Entity1_name': 'Gene'}, inplace=True)
                pharmgkb_drugs = pharmgkb_drugs[(pharmgkb_drugs['Entity1_type'] == 'Gene') & (pharmgkb_drugs['Association'] == 'associated')]
                pharmgkb_drugs['Entity2_name'] = pharmgkb_drugs['Entity2_name'].str.lower()
                pharmgkb_disease_drugs = pharmgkb_drugs[(pharmgkb_drugs['Entity2_name'] == self.disease)]
                pharmgkb_disease_drugs = pharmgkb_disease_drugs[['Gene']]
                pharmgkb_disease_drugs['pharmgkb_drugged_gene'] = 1
                gwascatalog = pd.read_table('./exploromics/data/gwas_catalog_v1.0-associations_e108_r2022-12-21.tsv', usecols = ['REPORTED GENE(S)', 'DISEASE/TRAIT', 'P-VALUE'], sep='\t')
                gene_chars = pd.read_table('./exploromics/data/hg19Rel92_AllgeneTypes_0kb.txt', sep='\t', header=None)
                gene_chars.rename(columns={0: 'Chr', 1:'Start', 2:'End', 3:'Gene', 4:'Gene_type', 5:'Source',}, inplace=True)
                gene_chars['Gene_length'] = gene_chars['End'] - gene_chars['Start']
                gene_chars = gene_chars[['Gene', 'Chr', 'Gene_length', 'Gene_type']]
                gene_chars['Chr'] = gene_chars['Chr'].astype('Int64')
                df_gene_chars = pd.merge(self.background_genes, gene_chars, on='Gene', how='left')
                gene_length_tmp = df_gene_chars[["Gene", "Chr", "Gene_length", "Gene_type"]]
                gene_length = gene_length_tmp.drop_duplicates(subset=['Gene'])

                if self.disease is not None and self.disease_id is not None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, self.open_targets_genes, pharmgkb_disease_drugs,
                                                                                    dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis])
                elif self.disease is not None and self.disease_id is None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, pharmgkb_disease_drugs,
                                                                            dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis])
                elif self.disease is None and self.disease_id is not None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, self.open_targets_genes, dgidb_drugs,
                                                                            hipred, sdi, pli, o_glncnac, rvis])

                elif self.disease is None and self.disease_id is None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, dgidb_drugs,
                                                                            hipred, sdi, pli, o_glncnac, rvis])

                
                if self.features_file is not None and self.disease is not None and self.disease_id is not None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, self.open_targets_genes, pharmgkb_disease_drugs,
                                                                                    dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis, self.features_file])
                    
                elif self.features_file is not None and self.disease is None and self.disease_id is not None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, self.open_targets_genes, dgidb_drugs,
                                                                                hipred, sdi, pli, o_glncnac, rvis, self.features_file])

                    
                elif self.features_file is not None and self.disease is not None and self.disease_id is None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, pharmgkb_disease_drugs,
                                                                                    dgidb_drugs, hipred, sdi, pli, o_glncnac, rvis, self.features_file])

                    
                elif self.features_file is not None and self.disease is None and self.disease_id is None:
                    background_df_eda = reduce(lambda x,y: pd.merge(x,y, on='Gene', how='left'), [self.background_genes, gene_length, gtex, dgidb_drugs,
                                                                                hipred, sdi, pli, o_glncnac, rvis, self.features_file]) 

                    
                if self.features_file is not None:
                    background_initial_columns = set(background_df_eda.columns)
                    background_df_eda_deduplicated = background_df_eda.drop_duplicates(subset=background_df_eda.columns, keep=False)
                    background_dropped_columns = background_initial_columns - set(background_df_eda_deduplicated.columns)
                    if background_dropped_columns:
                        print("The following background gene annotations were dropped (input features match automated integrated features):", background_dropped_columns)
                    background_df_eda = background_df_eda.drop_duplicates(subset=background_df_eda.columns, keep=False)

                else:
                    pass
                
                background_df_eda = background_df_eda.drop_duplicates(subset=['Gene'])
                background_df_eda.to_csv('./exploromics_results/background_genes/Background_genes_multi-omic_integration_output.csv', index=False) 
                self.background_df_eda = background_df_eda 
                background_gene_length_dict = defaultdict(list)

                for col in background_df_eda.columns[4:]:
                    corr = background_df_eda[str(col)].corr(background_df_eda['Gene_length'], method='spearman')
                    background_gene_length_dict[str(col)].append(corr)
                    background_gene_length_corr = pd.DataFrame.from_dict(background_gene_length_dict)
                    background_gene_length_corrT = background_gene_length_corr.T
                    background_gene_length_corrT.rename(columns={0: 'Gene length pearson correlation'}, inplace=True)
                background_gene_length_corrT.to_csv('./exploromics_results/background_genes/Background_gene_length_correlation.csv', index=False)
                print("Multi-omic data integrated for background genes...")
                return self.background_df_eda
        else:
            self.background_df_eda = None
            return self.background_df_eda

    def run(self):
        self.opentargets_association()
        self.gene_annotations()
        self.background_annotations()
    

if __name__ == '__main__':
    annotate = Annotation()
    annotate.run()