The appendix of my PhD thesis contained an analysis comparing the discriminatory power of the SNP panel created in this PhD thesis versus a panel of 15 microsatellite loci for the monitoring of *Vespa velutina* in the UK.
In order to compare the discriminatory power, I used a dataset of the same 86 individuals from across 2016-2024 for both SNPs and Microsatellites. I used the genotype data for both sets of genetic markers and did the below analysis.

## Related R package
* I used the related R package to compare the discriminatory power of SNPs vs Microsatellites for determining kin relationships.
* The genotype data for the real datasets was inputting into R and then was used to assign pair-wise relatedness coefficients to each pair of individuals.
* This was then used to simulate a dataset of 1,000 individuals per kin category (unrelated, half-sibling, full-sibling and parent-offspring) for both genetic markers.
* The dataset was then plotted on a histogram style graph to compare the discriminatory power.
* The code I used for the SNPs can be found under ```Related_SNPs.R``` and the code for microsatellites can be found under ```Related_Microsatellites.R```

## GenAlEx Excel Plug In
* 
