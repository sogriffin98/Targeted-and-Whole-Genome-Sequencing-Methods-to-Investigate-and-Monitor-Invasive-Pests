# SNP Selection
## Prepare Sample Reads and Compare to Reference Genome
### Trim sample Illumina reads using trim galore on --paired setting:
```
trim_galore --paired ALD_L4_S2_R1_001.fastq.gz ALD_L4_S2_R2_001.fastq.gz
trim_galore --paired AS_P8_S9_R1_001.fastq.gz AS_P8_S9_R2_001.fastq.gz
trim_galore --paired D8_P7_S1_R1_001.fastq.gz D8_P7_S1_R2_001.fastq.gz
trim_galore --paired GOS_L2_S4_R1_001.fastq.gz GOS_L2_S4_R2_001.fastq.gz
trim_galore --paired Jersey_rep_S6_R1_001.fastq.gz Jersey_rep_S6_R2_001.fastq.gz
trim_galore --paired TET_P6_S3_R1_001.fastq.gz TET_P6_S3_R2_001.fastq.gz
trim_galore --paired WOOL_P8_S7_R1_001.fastq.gz WOOL_P8_S7_R2_001.fastq.gz
```
### Index reference genome using samtools and bwa:
```
samtools faidx GCF_912470025.1_iVesVel2.1_genomic.fna
bwa index GCF_912470025.1_iVesVel2.1_genomic.fna.fai
```
### Align the reads:
```
cd ~/day1practical2/
mkdir -p ./bwa/sam/ 
mkdir -p ./bwa/bam/
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/ALD_L4_S2_R1_001_val_1.fastq.gz ./fastq/ALD_L4_S2_R2_001_val_2.fastq.gz > ./bwa/sam/ALD.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/AS_P8_S9_R1_001_val_1.fastq.gz ./fastq/AS_P8_S9_R2_001_val_2.fastq.gz > ./bwa/sam/AS.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/D8_P7_S1_R1_001_val_1.fastq.gz ./fastq/D8_P7_S1_R2_001_val_2.fastq.gz > ./bwa/sam/D8.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/GOS_L2_S4_R1_001_val_1.fastq.gz ./fastq/GOS_L2_S4_R2_001_val_2.fastq.gz > ./bwa/sam/GOS.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/Jersey_rep_S6_R1_001_val_1.fastq.gz ./fastq/Jersey_rep_S6_R2_001_val_2.fastq.gz > ./bwa/sam/JerseyRep.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/TET_P6_S3_R1_001_val_1.fastq.gz ./fastq/TET_P6_S3_R2_001_val_2.fastq.gz > ./bwa/sam/TET.sam
bwa mem ./GCF_912470025.1_iVesVel2.1_genomic.fna ./fastq/WOOL_P8_S7_R1_001_val_1.fastq.gz ./fastq/WOOL_P8_S7_R2_001_val_2.fastq.gz > ./bwa/sam/WOOL.sam
```
### Convert sam file into bam file:
```
samtools view -bh ./bwa/sam/ALD.sam > ./bwa/bam/ALD.bam
samtools view -bh ./bwa/sam/AS.sam > ./bwa/bam/AS.bam
samtools view -bh ./bwa/sam/D8.sam > ./bwa/bam/D8.bam
samtools view -bh ./bwa/sam/GOS.sam > ./bwa/bam/GOS.bam
samtools view -bh ./bwa/sam/JerseyRep.sam > ./bwa/bam/JerseyRep.bam
samtools view -bh ./bwa/sam/TET.sam > ./bwa/bam/TET.bam
samtools view -bh ./bwa/sam/WOOL.sam > ./bwa/bam/PORT.bam
```
### Sort the bam to make an index:
```
cd ~/heterozygosity/bwa/bam/
samtools sort -o ALD_sorted.bam ALD.bam
samtools index ALD_sorted.bam
samtools sort -o AS_sorted.bam AS.bam
samtools index AS_sorted.bam
samtools sort -o D8_sorted.bam D8.bam
samtools index D8_sorted.bam
samtools sort -o GOS_sorted.bam GOS.bam
samtools index GOS_sorted.bam
samtools sort -o JerseyRep_sorted.bam JerseyRep.bam
samtools index JerseyRep_sorted.bam
samtools sort -o TET_sorted.bam TET.bam
samtools index TET_sorted.bam
samtools sort -o WOOL_sorted.bam WOOL.bam
samtools index WOOL_sorted.bam
```
### Filter the alignment:
```
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in ALD_sorted.bam -out ALD_sorted.mq30.maxinsert.primaryalignment.bam
samtools index ALD_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in AS_sorted.bam -out AS_sorted.mq30.maxinsert.primaryalignment.bam
samtools index AS_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in D8_sorted.bam -out D8_sorted.mq30.maxinsert.primaryalignment.bam
samtools index D8_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in GOS_sorted.bam -out GOS_sorted.mq30.maxinsert.primaryalignment.bam
samtools index GOS_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in JerseyRep_sorted.bam -out JerseyRep_sorted.mq30.maxinsert.primaryalignment.bam
samtools index JerseyRep_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in TET_sorted.bam -out TET_sorted.mq30.maxinsert.primaryalignment.bam
samtools index TET_sorted.mq30.maxinsert.primaryalignment.bam
bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in WOOL_sorted.bam -out WOOL_sorted.mq30.maxinsert.primaryalignment.bam
samtools index WOOL_sorted.mq30.maxinsert.primaryalignment.bam
```
### Mark PCR duplicates using Picard:
```
picard MarkDuplicates I=ALD_sorted.mq30.maxinsert.primaryalignment.bam O=ALD_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=ALD_marked_dup_metrics.txt
picard MarkDuplicates I=AS_sorted.mq30.maxinsert.primaryalignment.bam O=AS_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=AS_marked_dup_metrics.txt
picard MarkDuplicates I=D8_sorted.mq30.maxinsert.primaryalignment.bam O=D8_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=D8_marked_dup_metrics.txt
picard MarkDuplicates I=GOS_sorted.mq30.maxinsert.primaryalignment.bam O=GOS_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=GOS_marked_dup_metrics.txt
picard MarkDuplicates I=JerseyRep_sorted.mq30.maxinsert.primaryalignment.bam O=JerseyRep_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=JerseyRep_marked_dup_metrics.txt
picard MarkDuplicates I=TET_sorted.mq30.maxinsert.primaryalignment.bam O=TET_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=TET_marked_dup_metrics.txt
picard MarkDuplicates I=WOOL_sorted.mq30.maxinsert.primaryalignment.bam O=WOOL_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam M=WOOL_marked_dup_metrics.txt
```
### Check quality of alignment:
```
mkdir qualimap
qualimap bamqc --java-mem-size=8g -bam ALD_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile ALD.pdf
qualimap bamqc --java-mem-size=8g -bam AS_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile AS.pdf
qualimap bamqc --java-mem-size=8g -bam D8_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile D8.pdf
qualimap bamqc --java-mem-size=8g -bam GOS_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile GOS.pdf
qualimap bamqc --java-mem-size=8g -bam JerseyRep_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile JerseyRep.pdf
qualimap bamqc --java-mem-size=8g -bam TET_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile TET.pdf
qualimap bamqc --java-mem-size=8g -bam WOOL_sorted.mq30.maxinsert.primaryalignment_marked_duplicates.bam -outdir ./qualimap/asianhornet_qualimap -outfile WOOL.pdf
```
### Variant detection:
```
mkdir -p ./vcf
bcftools mpileup -f ./asianhornetgenome/GCF_912470025.1_iVesVel2.1_genomic.fna --annotate AD,DP ./bwa/bam/*primaryalignment_marked_duplicates.bam | bcftools call -m -v -f GQ --skip-variants indels -o ./vcf/AH_subset_7inds.vcf
```
### Check the number of loci and individuals using vcftools:
```
vcftools --vcf AH_subset_7inds.vcf
```
### Filter by genotype:
```
vcftools --vcf AH_subset_7inds.vcf --minDP 6 --minGQ 18 --recode --out AH_subset_7inds_minDP6GQ18
```
## SNP Filtering:
### Biallelic Loci:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18.recode.vcf --max-alleles 2 --min-alleles 2 --recode --out AH_subset_7inds_minDP6GQ18_biallelic
```
### Missing Data:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic.recode.vcf --max-missing 0.14 --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing
```
### Minor Allele Frequency/Count:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing.recode.vcf --mac 2 --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2
```
### Maximum Depth Mean:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2.recode.vcf --site-mean-depth --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2
```
### Open file in text editor:
```
head AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2.ldepth.mean
```
### Use R to calculate mean and SD:
```
R
ldepth<-read.table("AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2.ldepth.mean", sep="\t", header = TRUE)
mean(ldepth$MEAN_DEPTH)
sd(ldepth$MEAN_DEPTH)
exit()
```
### Apply SD to VCF file (make sure to edit max-mean DP based on Mean):
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2.recode.vcf --max-meanDP 44 --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missinG_mac2_maxDP
```
### Check Coverage:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP.recode.vcf --depth --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP
```
### Filter for no missing:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP.recode.vcf --max-missing 1 --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing
```
### Linkage Disequilibrium (Thin 10k)
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing.recode.vcf --thin 10000 --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k
```
### Filter by Frequency
First I filtered by 0.5 frequency but this was too stringent:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k.recode.vcf --positions 0.5_freq_loci.txt --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.5freq
```
so I changed the filter to 0.4 - 0.6:
```
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k.recode.vcf --positions 0.40.6loci.txt --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq
```
### Filter by independant pairwise:
```
plink --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero.recode.vcf --indep-pairwise 10 5 0.5 --double-id --allow-extra-chr --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned
```
### Linkage Disequilibrium:
```
awk 'BEGIN{OFS="\t"} !/#/ {sub(/\./, $1"_"$2, $3)}1' AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero.recode.vcf > AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero.anotated.vcf

plink --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero.anotated.vcf --indep-pairwise 10 5 0.5 --double-id --allow-extra-chr --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned

sed 's/_/\t/2' AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned.prune.in > Good_LD_Loci.txt

vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero.recode.vcf --positions Good_LD_Loci.txt --recode --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned
```
## Selecting SNPs
After all the filters you should be left with 1011 SNP loci. I generated a python script following instructions from: https://www.biostars.org/p/334253/ to extract the SNPs from the reference genome and 100 flanking bases either side. 
```python3 flankingsequences.py```
