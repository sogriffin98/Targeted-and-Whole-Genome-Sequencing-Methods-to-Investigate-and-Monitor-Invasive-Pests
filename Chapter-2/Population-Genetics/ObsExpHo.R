# Set working directory
setwd("~/PhD/Thesis/Corrections")

# Load packages
library(vcfR)
library(adegenet)

# Function to calculate Ho, He, and FIS
calc_het <- function(vcf_file) {
  
  # Read VCF
  vcf <- read.vcfR(vcf_file)
  
  # Convert to genind object
  genind_obj <- vcfR2genind(vcf)
  
  # ---- Heterozygosity ----
  Ho <- mean(summary(genind_obj)$Hobs, na.rm = TRUE)
  He <- mean(summary(genind_obj)$Hexp, na.rm = TRUE)
  
  # ---- Inbreeding coefficient ----
  FIS <- 1 - (Ho / He)
  
  return(c(Ho = Ho, He = He, FIS = FIS))
}

# ---- Your VCF files ----
vcf_files <- c(
  "AH_104inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf", # All
  "AH_89inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf", # UK
  "AH_15inds_minDP6GQ18_0.7missing_mac2_biallelic_thinned1kb.recode.vcf"   # Europe
)

dataset_names <- c("All", "UK", "Europe")

# ---- Run analysis ----
results_list <- lapply(vcf_files, calc_het)

# ---- Format results ----
results_df <- data.frame(
  Dataset = dataset_names,
  Ho  = sapply(results_list, function(x) x["Ho"]),
  He  = sapply(results_list, function(x) x["He"]),
  FIS = sapply(results_list, function(x) x["FIS"])
)

print(results_df)






