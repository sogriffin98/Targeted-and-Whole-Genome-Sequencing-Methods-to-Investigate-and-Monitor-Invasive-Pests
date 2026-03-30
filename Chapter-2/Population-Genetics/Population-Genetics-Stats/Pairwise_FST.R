# Set working directory
setwd("~/PhD/Thesis/Corrections")

# Load packages
library(stringr)
library(vcfR)
library(SNPRelate)
library(SeqArray)
library(dplyr)

# Load metadata
meta <- read.csv("metadata.csv")

# Read VCF and clean sample names
vcf <- read.vcfR("AH_104inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf")

# Extract original names
orig_names <- colnames(vcf@gt)[-1]  # skip FORMAT column

# Clean them
clean_names <- orig_names %>%
  str_replace("^\\.\\/bam_104inds\\/", "") %>%   # remove ./bam_104inds/
  str_replace("^\\.\\/bam\\/", "") %>%           # remove ./bam/ if present
  str_replace("_combined", "") %>%
  str_replace("\\.mq30\\.maxinsert\\.primaryalignment\\.bam$", "")

# Apply cleaned names
colnames(vcf@gt) <- c("FORMAT", clean_names)

setdiff(clean_names, meta$Sample)
setdiff(meta$Sample, clean_names)

# Load metadata and check for mismatches
meta <- read.csv("metadata.csv")
setdiff(clean_names, meta$Sample)
setdiff(meta$Sample, clean_names)

# Convert VCF to GDS
write.vcf(vcf, file="cleaned.vcf")
snpgdsVCF2GDS("cleaned.vcf", "data.gds")
genofile <- snpgdsOpen("data.gds")

# Create population vector using metadata
popmap <- meta %>%
  filter(Sample %in% clean_names) %>%
  arrange(match(Sample, clean_names))

pop_vector <- factor(popmap$Country)

length(pop_vector)
length(clean_names)

# Compute global FST (Weir and Cockerham)
fst_global <- snpgdsFst(
  genofile,
  population = pop_vector,
  autosome.only = FALSE,
  method = "W&C84"
)

fst_global$Fst
