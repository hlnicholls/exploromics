# Exploromics
Tool for genetic exploratory data analysis and multi-omic data integration for machine learning.

- Intergrates multi-omic data for an input gene list
- Visualises gene characteristics (chromosome distribution, gene types, gene length correlation)
- Investigates annotated features (missingness and correlation) 
- Integrates multi-omic data for a background gene list to compare with input gene list (Kolmogorov–Smirnov testing)

## Introduction:

Machine learning applied to genomic data can have several pitfalls when the potentially biasing biological characteristics underlying the data have not been considered[1]. 

```Exploromics``` aims to assess a gene-level multi-omic dataset, providing an analysis of the biological factors that may impact machine learning applications. 

## Installation:
Requirements: Python3 (tested with v3.9.16)

```
conda env create -f exploromics.yml
conda activate exploromics
python setup.py install
```


## Run

- Set working directory to package location, then run:

```
exploromics --genelist [hgnc gene list text file] --features [features text file] --disease [disease name] --diseaseID [disease id (as used by OpenTargets)] --background [background gene list text file]
```

Example run:

```
exploromics --genelist example/genes.txt --features example/features.txt --disease 'Essential hypertension' --diseaseID 'EFO_0000537' --background example/other_genes.txt

```

## Arguments

```
required:
  --genelist GENELIST             Path to the gene list text file

optional arguments:
  -h, --help                      show this help message and exit
  --features [FEATURES]           Path to user features text file
  --disease [DISEASE]             Disease name
  --diseaseID [DISEASEID]         Disease ID (as used by OpenTargets)
  --background [BACKGROUND GENES] Background gene list text file
  ```
  
## Reference
1. Whalen, S., Schreiber, J., Noble, W.S. et al. Navigating the pitfalls of applying machine learning in genomics. Nat Rev Genet 23, 169–181 (2022). https://doi.org/10.1038/s41576-021-00434-9
