# Exploromics
Tool for genetic exploratory data analysis and multi-omic data integration for machine learning.

- Intergrates multi-omic data for an input genevlist
- Visualises gene characteristics (chromosome distribution, gene types, gene length correlation)
- Investigates annotated features (missingness and correlation) 
- Integrates multi-omic data for a background gene list to compare with input gene list 

## Installation:
Requirements: Python3 (tested with v3.9.16)

```
conda env create -f exploromics.yml
conda activate exploromics
python setup.py install
```


## Run

- Set command-line working directory to package location, then run:

```
exploromics --genelist [hgnc gene list text file] --features [features text file] --disease [disease name] --diseaseID [disease id (as used by OpenTargets)] --background [background gene list text file]
```

Example run:

```
exploromics --genelist example/genes.txt --disease 'Essential hypertension' --diseaseID 'EFO_0000537' --background example/other_genes.txt

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