### ============================================================
###  PH.D. CORRECTIONS: HETEROZYGOSITY, FIS, HWE, GLOBAL FST
###  Author: Sarah
### ============================================================

setwd("~/PhD/Thesis/Corrections")

### -------------------------
### Load packages
### -------------------------
library(vcfR)
library(adegenet)
library(pegas)
library(SNPRelate)
library(SeqArray)
library(dplyr)
library(stringr)

### ============================================================
### 1. FUNCTION: Fix genind objects (critical for GT-seq SNPs)
### ============================================================
fix_genind <- function(genind_obj) {
  
  # Ensure locus names exist
  if (is.null(locNames(genind_obj))) {
    locNames(genind_obj) <- paste0("Locus_", seq_len(nLoc(genind_obj)))
  }
  
  # Ensure allele names exist (simple placeholder if missing)
  if (is.null(alleles(genind_obj))) {
    alleles(genind_obj) <- rep(list(c("A","B")), nLoc(genind_obj))
  }
  
  # Remove loci with no genotype calls
  good_loci <- which(colSums(!is.na(tab(genind_obj))) > 0)
  genind_obj <- genind_obj[, good_loci]
  
  return(genind_obj)
}

### ============================================================
### 2. FUNCTION: Calculate Ho, He, FIS
### ============================================================
calc_het <- function(vcf_file) {
  vcf <- read.vcfR(vcf_file)
  genind_obj <- vcfR2genind(vcf)
  genind_obj <- fix_genind(genind_obj)
  
  s <- summary(genind_obj)
  Ho <- mean(s$Hobs, na.rm = TRUE)
  He <- mean(s$Hexp, na.rm = TRUE)
  FIS <- 1 - (Ho / He)
  
  return(list(genind = genind_obj, Ho = Ho, He = He, FIS = FIS))
}

### ============================================================
### 3. LOAD DATASETS
### ============================================================
vcf_files <- c(
  "AH_104inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf", # All
  "AH_89inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf", # UK
  "AH_15inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf"  # Europe
)

dataset_names <- c("All", "UK", "Europe")

results <- lapply(vcf_files, calc_het)
names(results) <- dataset_names

### ============================================================
### 4. ASSIGN POPULATIONS FROM METADATA (for FST / other uses)
### ============================================================
meta <- read.csv("metadata.csv")

assign_pops <- function(genind_obj) {
  inds <- indNames(genind_obj)
  submeta <- meta[match(inds, meta$Sample), ]
  genind_obj@pop <- factor(submeta$Country)  # "UK" / "Europe"
  return(genind_obj)
}

results$All$genind    <- assign_pops(results$All$genind)
results$UK$genind     <- assign_pops(results$UK$genind)
results$Europe$genind <- assign_pops(results$Europe$genind)

### ============================================================
### 5. SUMMARY: Ho, He, FIS
### ============================================================
het_df <- data.frame(
  Dataset = dataset_names,
  Ho = sapply(results, function(x) x$Ho),
  He = sapply(results, function(x) x$He),
  FIS = sapply(results, function(x) x$FIS)
)

print(het_df)

### ============================================================
### 6. HARDY–WEINBERG TESTS (chi-square, one population per dataset)
### ============================================================
run_hwe <- function(genind_obj) {
  # Treat each dataset as a single panmictic population
  pop(genind_obj) <- factor("OnePop")
  hwe <- hw.test(genind_obj, B = 0)   # chi-square test
  pvals <- as.numeric(hwe)
  sig <- sum(pvals < 0.05, na.rm = TRUE)
  total <- length(pvals)
  return(list(pvals = pvals, sig = sig, total = total))
}

hwe_results <- lapply(results, function(x) run_hwe(x$genind))

for (i in 1:length(dataset_names)) {
  cat("\n---", dataset_names[i], "---\n")
  cat("Significant loci:", hwe_results[[i]]$sig, "out of", hwe_results[[i]]$total, "\n")
  cat("Proportion:", hwe_results[[i]]$sig / hwe_results[[i]]$total, "\n")
}

### ============================================================
### 7. GLOBAL FST (Weir & Cockerham) + PERMUTATION TEST
### ============================================================

### Clean sample names in ALL dataset VCF
vcf_all <- read.vcfR(vcf_files[1])
orig_names <- colnames(vcf_all@gt)[-1]

clean_names <- orig_names %>%
  str_replace("^\\.\\/bam_104inds\\/", "") %>%
  str_replace("^\\.\\/bam\\/", "") %>%
  str_replace("_combined", "") %>%
  str_replace("\\.mq30\\.maxinsert\\.primaryalignment\\.bam$", "")

colnames(vcf_all@gt) <- c("FORMAT", clean_names)

### Match metadata order
popmap <- meta %>%
  filter(Sample %in% clean_names) %>%
  arrange(match(Sample, clean_names))

pop_vector <- factor(popmap$Country)  # "UK" / "Europe"

### Convert to GDS
write.vcf(vcf_all, file = "cleaned.vcf")
snpgdsVCF2GDS("cleaned.vcf", "data.gds")
genofile <- snpgdsOpen("data.gds")

### Compute global FST
fst_global <- snpgdsFst(
  genofile,
  population = pop_vector,
  autosome.only = FALSE,
  method = "W&C84"
)

cat("\nGlobal FST:", fst_global$Fst, "\n")

### Permutation test for FST (shuffle UK/Europe labels)
perm_fst <- replicate(
  1000,
  {
    shuffled <- sample(pop_vector)
    snpgdsFst(genofile, population = shuffled, autosome.only = FALSE)$Fst
  }
)

p_value <- mean(perm_fst >= fst_global$Fst)

cat("Permutation P-value:", p_value, "\n")
