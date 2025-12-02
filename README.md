# Genomic Tools to Enable Epidemiological Monitoring of Eukaryotic Pests
This GitHub repository contains all of the scripts and code used during my PhD project titled "Genomic Tools to Enable Epidemiological Monitoring of Eukaryotic Pests" at Newcastle University from September 2022 to July 2026. The code is split into the relevant chapters in the folders of this repository to allow ease of reproducibility. The study species for this research were *Vespa velutina* (yellow legged hornet) and *Meloidogyne fallax* (false Columbia root knot nematode). *V. velutina* is an invasive insect species in the UK which poses a threat to apiculture. *M. fallax* is a plant pest which burrows into the roots of crops causing external galling and necrosis. 

## Chapter 1: Introduction Chapter
* In the introduction, I created a map of the UK with the individual insects and nests found for *Vespa velutina*. 
* I used R studio to create this map based on location data from Defra: https://www.gov.uk/government/publications/asian-hornet-uk-sightings/asian-hornet-sightings-recorded-since-2016. 
* The script and Microsoft Excel file is in the Literature-Review folder of this repository: ```Mapping_of_Hornets.R``` and ```mapping_of_hornets.xlsx```

## Chapter 2: Development of a Genotyping in Thousands by Sequencing (GT-Seq) Assay and Population Genetics Analysis of *Vespa velutina* (Yellow Legged Hornet) 
The first data chapter of my PhD project focuses on *Vespa velutina* (yellow legged or Asian hornet). The first part of the chapter discusses the process involved with the selection of SNP loci and the primer design. Then the second part checks the data entegrity of the Illumina MiSeq data before it is used for data analysis. The data analysis part of this chapter is split into the 2 parts, a and b, below (Population Genetics and Kinship Analysis) where I have provided the input files and code used for each analysis type in separate folders.
### SNP Selection and Primer Design
#### 1. SNP Selection
This folder contains python scripts and commandline code used for:
* Trimming WGS data
* Aligning WGS to the Reference Genome
* Filtering the SNPs
#### 2. Primer Design
This folder contains the python and perl scripts and commandline code used for:
* N-masking the genome
* Primer Design using Primer3 commandline
* Checking the primers using the GT-seq primer check script
* Filtering of primers by GC% and melting temperature
### Illumina MiSeq Data Entegrity


NEED TO SORT!!!!!!!!!!!!!!!!!!!!!


### Data Chapter 1a: Population Genetics Analysis
#### 1. Principal-Components-Analysis
This folder contains R code and input files for:
* UK 1 per nest analysis: UK_1pernest_PCA.R, UK_pca_results.eigenval, UK_pca_results.eigenvec, sample_ID_UK_1pernest.csv
* UK/EU 1 per nest analysis: UKEU_1pernest_PCA.R, UKEU_pca_results.eigenval, UKEU_pca_results.eigenvec, sample_ID_UKEU_1pernest.csv
#### 2. STRUCTURE
This folder contains R code and input files for:
* UK 1 per nest analysis: UK_1pernest_STRUCTURE.R, structure_summary_uk.csv
* UK/EU 1 per nest analysis: UKEU_1pernest_STRUCTURE.R, structure_summary_uk_eu.csv
#### 3. Isolation-By-Distance
This folder contains R code and input files for:
* UK 1 per nest analysis: UK_1pernest_IBD.R, UK_1pernest_coords.csv, UK_1pernest_data_plink_output.raw
* UK/EU 1 per nest analysis: UKEU_1pernest_IBD.R, UKEU_1pernest_coords.csv, UKEU_1pernest_data_plink_output.raw
### Data Chapter 1b: Kinship Analysis
#### 1. Related

## Data Chapter 2: Genome assembly of *Meloidogyne fallax*
* In this data chapter, I created a draft reference genome for *Meloidogyne fallax* or false Columbia root knot nematode.
* The genome was published in the Journal of Nematology: https://doi.org/10.2478/jofnem-2025-0016
* The associated scripts and code is in the Data-Chapter-2 folder of this repository.
* The reference genome can be found here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048569365.1/
* The raw Illumina and Oxford Nanopore Technologies NGS data can be found here: https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=44309162

## Appendix: Comparison between SNPs vs Microsatellites for 
* In the appendix, I did a comparison of the discrimatory power of SNPs vs Microsatellites for *Vespa velutina* (yellow legged hornet).
* Microsatellites are currently used by Fera Science Limited to monitor yellow legged hornets in the UK
* I compared the two markers using the related R package and GenAlEx to better understand the ability to distinguish between kin groups
