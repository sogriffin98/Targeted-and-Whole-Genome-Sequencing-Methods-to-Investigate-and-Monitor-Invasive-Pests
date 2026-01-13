The appendix of my PhD thesis contained an analysis comparing the discriminatory power of the SNP panel created in this PhD thesis versus a panel of 15 microsatellite loci for the monitoring of *Vespa velutina* in the UK.
In order to compare the discriminatory power, I used a dataset of the same 86 individuals from across 2016-2024 for both SNPs and Microsatellites. I used the genotype data for both sets of genetic markers and did the below analysis.

## Related R package
* I used the related R package to compare the discriminatory power of SNPs vs Microsatellites for determining kin relationships.
* The genotype data for the real datasets was inputting into R and then was used to assign pair-wise relatedness coefficients to each pair of individuals.
* This was then used to simulate a dataset of 1,000 individuals per kin category (unrelated, half-sibling, full-sibling and parent-offspring) for both genetic markers.
* The dataset was then plotted on a histogram style graph to compare the discriminatory power.
* The code I used for the SNPs can be found under ```Related_SNPs.R``` and the code for microsatellites can be found under ```Related_Microsatellites.R```
* The raw data files for both can be found under: ```snps_raw_data.txt``` and ```Related_microsatellites_raw_data.txt```

## GenAlEx Excel Plug In
* I used GenAlEx (https://biology-assets.anu.edu.au/GenAlEx/Download.html) to compare the probability of identity and probability of identity of siblings using the SNP and Microsatellite data.
* Using the genotype data from the SNPs and Microsatellites I compared the two markers using the genotype data for each. The microsatellite data was previously completed by Fera Science and was copied into a new Microsoft Excel document for the 89 individuals. The SNP genotype data was copied from the STRUCTURE file for the 89 individuals
* Then in Excel click the GenAlEx tab > Frequency-Based > Multilocus > Prob. Identity
* Fill out the Parameters that match the data and click OK
* There should be a new sheet created with the probability of identity and probability of identity of siblings for each marker type
