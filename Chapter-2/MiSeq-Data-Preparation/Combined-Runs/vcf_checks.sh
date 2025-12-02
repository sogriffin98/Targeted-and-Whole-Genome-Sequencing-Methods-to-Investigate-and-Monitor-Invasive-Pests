#!/bin/bash

# Replace <> with your specific identifier
PREFIX="AH_400_inds"

# Apply filters for minimum depth and genotype quality
vcftools --vcf ${PREFIX}.vcf --minDP 6 --minGQ 18 --recode --out ${PREFIX}_minDP6GQ18

# Filter for biallelic variants
vcftools --vcf ${PREFIX}_minDP6GQ18.recode.vcf --max-alleles 2 --min-alleles 2 --recode --out ${PREFIX}_minDP6GQ18_biallelic

# Generate depth statistics
vcftools --vcf ${PREFIX}.vcf --depth --out ${PREFIX}_depth

# Calculate site mean depth
vcftools --vcf ${PREFIX}.vcf --site-mean-depth --out ${PREFIX}_meandepthpersite

# Calculate site quality
vcftools --vcf ${PREFIX}.vcf --site-quality --out ${PREFIX}_sitequality

# Filter based on missing data threshold
vcftools --vcf ${PREFIX}.vcf --max-missing 0.5 --recode --out ${PREFIX}_0.5missing

# Apply depth and missing data filters
vcftools --vcf ${PREFIX}.vcf --minDP 6 --max-missing 0.6 --recode --out ${PREFIX}_minDP6_0.6missing

# Generate missingness report
vcftools --vcf ${PREFIX}_minDP6_0.6missing.recode.vcf --missing-indv --out ${PREFIX}_0.6missing_mindepth6_MissingnessReport

# Calculate site mean depth for filtered data
vcftools --vcf ${PREFIX}_0.5missing.recode.vcf --site-mean-depth --out ${PREFIX}_0.5missing_MEANSITEDEPTH

# Print completion message
echo "VCF processing completed."
