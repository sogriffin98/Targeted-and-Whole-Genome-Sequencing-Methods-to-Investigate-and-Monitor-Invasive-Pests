# GT-seq_Hornets_MiSeq_Analysis
This is the data analysis I performed separately on the first few Illumina MiSeq runs to check that the SNP loci has been amplified and to remove or dilute any overrepresented SNPs from the primer pool. 

## FastQC
* First run FastQC analysis on each of the reads:
```
conda activate fastqc
fastqc *.fastq.gz
```
## MiSeq Run 1 Primer Drop Data Analysis Part 1 (27th January 2025)
I followed the instructions from Nate Campbell's GitHub page (https://github.com/GTseq/GTseek_utils/blob/Main/README_GTseek-MultiPCR-Analysis.txt) from line 15.
### 1. Unzip all the files 
```
gunzip -k *.fastq.gz
```
### 2. Cat all the forward and reverse together but leave the undetermined forward and reverse out
```
cat *L001_R1_001.fastq > R1.fastq
cat *L001_R2_001.fastq > R2.fastq
```
### 3. Combine the forward and reverse reads together using perl script
```
perl paste_fastq.pl R1.fastq R2.fastq > pasted.fq
```
### 4. Create a hash file from the combined read file
```
perl HashSeqs.pl pasted.fq HORNETS > file.hash
```
### 5. Then run the primer interaction test on that hash file. The 471_Primer_Seqs.txt file is a tab deliminated file containing Locus Name, forward sequence, reverse sequence
```
 perl GTseq_Primer-Interaction-Test_v3.pl 471_Primer_Seqs.txt file.hash > file_PI-test.txt
```
### 6. The output of the above command is a very large file with a TON of information, most of which is useless, so I suggest truncating the file like this:
```
grep -A20 '[Pp]rimer\|[Bb]lack' file_PI-test.txt > file_PI-test-trunc.txt
```
### 7. Repeat Step 3-6 for the undetermined forward and reverse
### 8. Check the truncated output for which primers to drop

## MiSeq Run 1 Primer Drop Data Analysis Part 2 (28th January 2025)
These steps were to check the outputted sequencing file for overrepresented primers
### 1. Use AdapterRemoval (https://github.com/MikkelSchubert/adapterremoval?tab=readme-ov-file) to trim and merge reads to create 7 files with all the individuals merged
```
conda activate adapterremoval
AdapterRemoval --file1 *_R1_001.fastq --file2 *_R2_001.fastq --trimns --trimqualities --collapse
```
### 2. Run fastqc
```
conda activate fastqc
fastqc your_output.collapsed
fastqc your_output.pair1.truncated
fastqc your_output.pair2.truncated
```
### 3. Open the FastQC data file in the zipped folder and copy across the overrepresented sequences to an excel sheet
### 4. Take the first 15 bases of each primer and compare them to the list of 480 SNP primers to look for potential primers which could be removed
### 5. Looking at percentages and discard primers which are above 0.5%

## Alignment to reference genome (30th January 2025)
### 1. Run AdapterRemoval on each sample individually
```
conda activate adapterremoval
./adapterremoval.sh
```
### 2. Rename the collapsed files in the AdapterRemoval folder
```
chmod +x rename.sh
./rename.sh
```
### 3. Index the reference genome (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_912470025.1/) and align each of the individuals
```
mkdir -p ./bwa/sam/ 
mkdir -p ./bwa/bam/
conda activate bwa
chmod +x bwa_samse_alignment.sh
./bwa_samse_alignment.sh
```
### 4. Convert sam to bam
```
conda activate samtools
chmod +x sam_to_bam.sh
./sam_to_bam.sh 
```
### 5. Sort and Index Bam Files
```
chmod +x bam_sort_index.sh
./bam_sort_index.sh 
```
### 6. Filter the alignment by mapping quality and primer alignment 
```
conda activate bamtools
chmod +x filter_bam.sh
./filter_bam.sh 
```
### 7. Mpileup
```
conda activate bcftools
bcftools mpileup -f ./reference/GCF_912470025.1_iVesVel2.1_genomic.fna --annotate AD,DP ./bam/*primaryalignment.bam | bcftools call -m -v -f GQ --skip-variants indels -o ./vcf/AH_48_inds.vcf
```
### 8. Checking vcf file
* I created a bash script for checking and filtering the vcf file:
```
chmod +x vcf_checks.sh
./vcf_checks.sh
```
* The above bash script uses these commands from this link: https://speciationgenomics.github.io/filtering_vcfs/
```
vcftools --vcf AH_<>_indivs.vcf --minDP 6 --minGQ 18 --recode --out AH_<>_indivs_minDP6GQ18.vcf
vcftools --vcf AH_<>_indivs_minDP6GQ18.recode.vcf --max-alleles 2 --min-alleles 2 --recode --out AH_<>_indivs_minDP6GQ18_biallelic
vcftools --vcf AH_<>_indivs.vcf --depth --out AH_<>_indivs_depth
vcftools --vcf AH_<>_indivs.vcf --site-mean-depth --out AH_<>_indivs_meandepthpersite
vcftools --vcf AH_<>_indivs.vcf --site-quality --out AH_<>_indivs_sitequality
```
* Filter for 50% missing data (also included in bash script):
```
vcftools --vcf AH_<>_inds.vcf --max-missing 0.5 --recode --out AH_<>_inds_0.5missing
```
* Filter for a minimum depth of 6 and 40% missing data (also included in bash script):
```
vcftools --vcf AH_<>_inds.vcf --minDP 6 --max-missing 0.6 --recode --out AH_<>_inds_minDP6_0.6missing
```
### 10. Open the VCF Files and analyse in Microsoft Excel to compare presence of SNP primers between original, 50% filtered and min depth 40% filter VCF files
### 11. Run below code to check which samples have coverage below 20%. Open the report, order the individuals by the final column and then repeat sequencing for individuals below 0.02
```
vcftools --vcf AH_<>_indivs_minDP6_0.6missing.recode.vcf --missing-indv --out AH_<>_indivs_0.5missing_mindepth6_MissingnessReport
```
### 12. Also run vcftools to calculate mean site depth for each loci to see which primers to dilute. Order final column and then dilute primers above 1000 x coverage:
```
vcftools --vcf AH_<>_inds_0.5missing.recode.vcf --site-mean-depth --out AH_<>_inds_0.5missing_MEANSITEDEPTH
```
