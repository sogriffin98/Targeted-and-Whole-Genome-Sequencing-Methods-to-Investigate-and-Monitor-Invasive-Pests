# Primer Design
## Convert .vcf file into .bed file
First both the original .vcf file needs to be converted into a .bed file:
```
conda activate vcftools
vcftools --vcf AH_subset_7inds.vcf --site-depth --out AH_subset_7inds_sitedepth 
vcftools --vcf AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned.recode.vcf --site-depth --out AH_subset_7inds_minDP6GQ18_biallelic_0.14missing_mac2_maxDP_nomissing_thin10k_0.4_0.6freq_exhetero_LDpruned_recode
```
Then open the file in Excel. Copy the second column and add 1 to each number. Copy the first 3 columns into a Notepad file and save the file with a .bed extension. This will be the .bed files needed for step 1 and 2 below.

## Whole genome data n-masking:
### Needed:
1. Fasta file with the whole genome
2. Bed file with positions of all the variants found in the genome
3. Bed file with name of all the loci of interest, position -101bp, position +100bp
4. Txt file of the alleles for each variants of interest formatted as: Loci ID, [Allele 1/Allele 2]
   
### Step 1: n-mask whole genome file using bedtools - this allows to have a copy of the masked genome already made, and to only have to repeat the following steps if more variants need to be extracted later on down the line
```
bedtools maskfasta -fi <genome.fa> -bed <allvariants.bed> -fo <masked_genome.fa>
```
### Step 2: extract 200 base pairs fragments with the variants of interest in the middle
```
bedtools getfasta -fi <masked_genome.fa> -bed <variants_of_interest.bed> -name -fo <masked_variants_of_interest.fa>
```
### Step 3: replace masked allele in each sequence using py_replace_alleles.py script. The output is a fasta file with nmasked variants.
```
python3 py_replace_alleles.py <masked_variants_of_interest.fa> <alleles_info.txt>
```
## Primer Design Using Primer3
### Fasta2primer3
Python script to output primers forward and backward from a fasta file using Primer3 (https://github.com/frba/fasta2primer3/blob/master/README.md)
```
conda deactivate
python fasta2primer3.py Primer3_input.fa
```
This python script doesn't fully work but it filters the invalid sequences out. The script outputs 2 FASTA files: invalidseqs and validseqs. The validseqs FASTA file can be converted into a BoulderIO file:

### Converting the fasta file into BoulderIO
Primer3 doesn't read fasta files as an input file. The tool requires BoulderIO files as an input. To make a BoulderIO file from a fasta file you use a programme called fastaq (https://github.com/sanger-pathogens/Fastaq/blob/master/README.md). I used the following commands:
```
fastaq to_boulderio Primer3_input.fa-validseqs.fasta Primer3_input-validseqs.io
```
### Primer3
Primer3 is a primer design programme which has both a web interface and commandline functionality (https://primer3.org/manual.html). In this case I used the commandline functionality as I was wanting to design primers for hundreds of SNP loci. 
On the web interace it has an option to download the settings. Download the settings file and save it as a text file called ```primer3_settings.txt```
Using the BoulderIO and settings text file I used the commandline primer3 capability to design the primers. I did this by:
```
primer3_core --format_output --default_version=2 --p3_settings_file=primer3_settings.txt --output=valid_primers.txt Primer3_input-validseqs.io
primer3_core --format_output --default_version=2 --p3_settings_file=primer3_settings.txt --output=valid_primers.csv Primer3_input-validseqs.io
primer3_core --default_version=2 --p3_settings_file=primer3_settings.txt --output=valid_primers2.csv Primer3_input-validseqs.io
```
### Preparing file for GT-seq Primer Check Step
Using Primer3 above you can selece .csv as a file output which is the file format needed for the GT-seq primer check step. However, the .csv file is not in the correct format as specified in the perl script. It requires the format to be: Name\tFWD-Primer\tREV-Primer in a .csv file.

So in order to change the format of the .csv to the required one you have to do the following:
* Open the valid_primer2.csv file
* Delete each row that is not the name and the first forward and reverse primer (i.e. the primer that is called 0)
* You should be left with a file with the name and the 2 primers
* Then this can be formatted into 3 tab deliminated columns: Name\tFWD-Primer\tREV-Primer

### GT-seq Multiplex PCR Primer Check Step
This step was done using a Perl script from the GT-seq GitHub (https://github.com/GTseq/GTseq-Pipeline).
```
GTseq_PrimerCheck.pl GTseq_PrimerCheck_Primers.csv
```
## Final Filtering
* The next filtering step applied was the melting temperature (Tm), the temperature at which half of the oligonucleotide is annealed to its compliment (Whitman, 2023).
* In the case of this study, the primers were filtered to ensure that they had a Tm of over 50 degrees and a Tm difference between left and right primer of less than 5 degrees (Premier Biosoft., 2023, Thermofisher Scientific., 2019).
* This filtering step did not remove any primers from the list.
* The final filtering step was to assess the GC content of each of primer.
* The GC content of primers should ideally be 40-60%. In the case of this study, primers with a GC% of between 50-61.11% were kept.
* This removed 1,212 primers leaving 482 primers.
* For ease of ordering 96-well primer plates, 2 more primers were removed to leave a final total of 480 primer pairs for use in the Genotyping in Thousands by Sequencing (GT-seq) assay. 

## Ordering the Primers
The primers were ordered from IDT in a 96 well format with the forward and reverse primer pairs within the same well to reduce the cost of the project.
