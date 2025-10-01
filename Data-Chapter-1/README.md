# Data Chapter 1
The first data chapter of my PhD project focuses on *Vespa velutina* (yellow legged or Asian hornet). The first part of the chapter discusses the process involved with the selection of SNP loci and the primer design. The data analysis part of this chapter is split into the 2 parts, a and b, below (Population Genetics and Kinship Analysis) where I have provided the input files and code used for each analysis type in separate folders.
## SNP Selection and Primer Design
### 1. SNP Selection
This folder contains python scripts and commandline code used for:
* Trimming WGS data
* Aligning WGS to the Reference Genome
* Filtering the SNPs
### 2. Primer Design
This folder contains the python and perl scripts and commandline code used for:
* N-masking the genome
* Primer Design using Primer3 commandline
* Checking the primers using the GT-seq primer check script
* Filtering of primers by GC% and melting temperature

## Data Chapter 1a: Population Genetics Analysis
### 1. Principal-Components-Analysis
This folder contains R code and input files for:
* **UK 1 per nest analysis**: UK_1pernest_PCA.R, UK_pca_results.eigenval, UK_pca_results.eigenvec, sample_ID_UK_1pernest.csv
* **UK/EU 1 per nest analysis**: UKEU_1pernest_PCA.R, UKEU_pca_results.eigenval, UKEU_pca_results.eigenvec, sample_ID_UKEU_1pernest.csv

### 2. STRUCTURE
This folder contains R code and input files for:
* **UK 1 per nest analysis**: UK_1pernest_STRUCTURE.R, structure_summary_uk.csv
* **UK/EU 1 per nest analysis**: UKEU_1pernest_STRUCTURE.R, structure_summary_uk_eu.csv

### 3. Isolation-By-Distance
This folder contains R code and input files for:
* **UK 1 per nest analysis**: UK_1pernest_IBD.R, UK_1pernest_coords.csv, UK_1pernest_data_plink_output.raw
* **UK/EU 1 per nest analysis**: UKEU_1pernest_IBD.R, UKEU_1pernest_coords.csv, UKEU_1pernest_data_plink_output.raw

## Data Chapter 1b: Kinship Analysis
### Related
