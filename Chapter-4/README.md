# *Meloidogyne fallax* Genome Assembly
Genome assembly steps used to assemble the Meloidogyne fallax genome using Illumina and Oxford Nanopore Technologies reads.

## Quality Control Raw Reads:
### FastQC:
#### Oxford Nanopore Technologies:
```
fastqc MF1.fastq.gz
```
#### Illumina:
```
fastqc MF1_S12_R1_001.fastq.gz
fastqc MF1_S12_R2_001.fastq.gz
```
### NanoStat:
#### Oxford Nanopore Technologies:
```
NanoStat --fastq MF1.fastq.gz -o /data/ssd2/sgriffin/fallax/ONT/QC/Raw/ -n MF1_StatReport
```
#### Illumina:
```
NanoStat --fastq MF1_S12_R1_001.fastq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Raw/ -n MF1_R1_StatReport
NanoStat --fastq MF1_S12_R2_001.fastq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Raw/ -n MF1_R2_StatReport
NanoStat --fastqc MF1.fastq.gz -o /data/ssd2/sgriffin/fallax/ONT/QC/Raw/ -n MF1_StatReport
```
## Trimming of raw reads:
### Trim Galore (Illumina):
```
trim_galore --paired MF1_S12_R1_001.fastq.gz MF1_S12_R2_001.fastq.gz
trim_galore --paired --hardtrim5 235 --hardtrim3 215 MF1_S12_R1_001.fastq.gz MF1_S12_R2_001.fastq.gz
```
### Prowler (Oxford Nanopore Technologies):
```
gunzip -d MF1.fastq.gz
python3 TrimmerLarge.py -f "MF1.fastq" -i "/data/ssd2/sgriffin/fallax/ONT/ONT_Raw_MF/” -o "/data/ssd2/sgriffin/fallax/ONT/" -w 100 -l 1000 -c LT -g U0 -m S -q 7 -d 0 -r .fastq
```
## Quality Control Trimmed Reads:
### FastQC:
#### Oxford Nanopore Technologies:
```
fastqc MF1TrimLT-U0-S7W100L1000R0.fastq
```
#### Illumina:
```
fastqc MF1_S12_R1_001_val_1.fastq.gz
fastqc MF1_S12_R2_001_val_2.fastq.gz
fastqc MF1_S12_R1_001.235bp_5prime.fq.gz
fastqc MF1_S12_R2_001.235bp_5prime.fq.gz
fastqc MF1_S12_R1_001.235bp_5prime.215bp_3prime.fq.gz
fastqc MF1_S12_R2_001.235bp_5prime.215bp_3prime.fq.gz
```
### NanoStat:
#### Oxford Nanopore Technologies:
```
NanoStat --fastq MF1TrimLT-U0-S7W100L1000R0.fastq -o /data/ssd2/sgriffin/fallax/ONT/QC/ -n MF1_Prowler_StatReport
```
#### Illumina:
```
NanoStat --fastq MF1_S12_R1_001_val_1.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R1_TrimGalore_StatReport
NanoStat --fastq MF1_S12_R2_001_val_2.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R2_TrimGalore_StatReport
NanoStat --fastq MF1_S12_R1_001.235bp_5prime.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R1_235bp5prime_StatReport
NanoStat --fastq MF1_S12_R2_001.235bp_5prime.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R2_235bp5prime_StatReport
NanoStat --fastq MF1_S12_R1_001.235bp_5prime.215bp_3prime.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R1_235bp5prime_215bp3prime_StatReport
NanoStat --fastq MF1_S12_R2_001.235bp_5prime.215bp_3prime.fq.gz -o /data/ssd2/sgriffin/fallax/Illumina/QC/Trim/ -n MF1_R2_235bp5prime_215bp3prime_StatReport
```

## Assembly of Genomes:
### Flye (Oxford Nanopore Technologies):
```
flye -o MF1_ONT_Flye --threads 16 --nano-raw MF1TrimLT-U0-S7W1000L1000R0.fastq -m 1000
```
### Spades (Illumina):
```
spades.py -1 MF1_S12_R1_001_trimmed.fq.gz -2 MF1_S12_R2_001_trimmed.fq.gz -o MF1_trimmed_Spades -t 40 -m 800
spades.py -1 MF1_S12_R1_001.235bp_5prime.fq.gz -2 MF1_S12_R2_001.235bp_5prime.fq.gz -o MF1_5prime235_Spades -t 40 -m 800
spades.py -1 MF1_S12_R1_001.235bp_5prime.215bp_3prime.fq.gz -2 MF1_S12_R2_001.235bp_5prime.215bp_3prime.fq.gz -o MF1_5prime235_3prime215_Spades -t 40 -m 800
```
### Hybrid 
```
spades.py -1 MF1_S12_R1_001_val_1.fq.gz -2 MF1_S12_R2_001_val_2.fq.gz --nanopore MF1_Prowler.fastq -o MF1_Prowler_TrimGaloreNormal_SPAdesHybrid -t 40 -m 800
spades.py -1 MF1_S12_R1_001.235bp_5prime.fq.gz -2 MF1_S12_R2_001.235bp_5prime.fq.gz --nanopore MF1_Prowler.fastq -o MF1_Prowler_TrimGalore5prime_SPAdesHybrid -t 40 -m 800
spades.py -1 MF1_S12_R1_001_val_1.fq.gz -2 MF1_S12_R2_001_val_2.fq.gz --nanopore MF1_Prowler.fastq -o MF1_Prowler_TrimGalore5prime3prime_SPAdesHybrid -t 40 -m 800
```
## Genome Assembly Quality Control:
### Busco:
#### Illumina:
```
busco -i contigs.fasta --auto-lineage -o TrimGalore_SPAdes_BUSCO -m genome
busco -i contigs.fasta --auto-lineage -o TrimGalore_SPAdes_BUSCO_5Prime_3Prime -m genome
busco -i contigs.fasta --auto-lineage -o TrimGalore_SPAdes_BUSCO_5Prime -m genome
mkdir BUSCO_summaries_Illumina/
cp OUT1/short_summary.*.lineage_odb10.OUT1.txt ./BUSCO_summaries/
cp OUT2/short_summary.*.lineage_odb10.OUT2.txt ./BUSCO_summaries/
generate_plot.py -wd BUSCO_summaries_Illumina
```
#### Oxford Nanopore Technologies:
```
busco -i /data/ssd2/sgriffin/fallax/ONT/MF1_ONT_Flye/assembly.fasta --auto-lineage -o MF1_ONT_Flye_Busco -m genome
mkdir BUSCO_summaries/
cp OUT1/short_summary.*.lineage_odb10.OUT1.txt ./BUSCO_Summaries_ONT/
cp OUT2/short_summary.*.lineage_odb10.OUT2.txt ./BUSCO_Summaries_ONT/
generate_plot.py -wd BUSCO_Summaries_ONT
```
#### Hybrid:
```
busco -i contigs.fasta --auto-lineage -o Hybrid_ProwlerNormal_BUSCO -m genome
busco -i contigs.fasta --auto-lineage -o Hybrid_Prowler5prime_BUSCO -m genome
busco -i contigs.fasta --auto-lineage -o Hybrid_Prowler5prime3prime_BUSCO -m genome
mkdir BUSCO_Summaries_Hybrid/
cp OUT1/short_summary.*.lineage_odb10.OUT1.txt ./BUSCO_Summaries_Hybrid/
cp OUT2/short_summary.*.lineage_odb10.OUT2.txt ./BUSCO_Summaries_Hybrid/
generate_plot.py -wd BUSCO_Summaries_Hybrid
```
### Assembly Statistics
```
assemblyStatistics -f contigs.fasta -l 100
assemblyStatistics contigs.fasta
```
or 
```
assembly-stats <>.fasta
```
### QUAST
```
quast.py <>.fasta -o <assemblyname>_quast
```
### bbmap/bbtools
```
conda activate bbmap
```
Statistics on the assembly
```
stats.sh in=contigs.fa
```
Statistics on multiple assemblies
```
statswrapper.sh in=a.fa,b.fa,c.fa format=6
```
Print GC and length information per sequence:
```
stats.sh in=contigs.fa gc=gc.txt gcformat=4
```

## Ragtag
Information about the tool can be found at the GitHub page:https://github.com/malonge/RagTag
### Correct
```
ragtag.py correct GCA_015183035.1_ASM1518303v1_genomic.fasta MF1_Illumina_SPAdes_TrimGalore_btkfiltered.fasta -R MF1_S12_R1_001_val_1.fq -R MF1_S12_R2_001_val_2.fq -T sr
```
### Scaffold
```
ragtag.py scaffold GCA_015183035.1_ASM1518303v1_genomic.fasta ragtag.correct.fasta
```
### Quast
```
quast.py ragtag.scaffold.fasta -o RagTag_patch_QUAST
```

## TGS Gap Closer
Information about the tool can be found at the GitHub page:https://github.com/BGI-Qingdao/TGS-GapCloser
```
tgsgapcloser --scaff ragtag.scaffold.fasta --reads MF1_ONT_Flye_assembly.fasta --output MF1_RagTagCorrScaff_TGSGapCloserONT --ne
```
### Quast
```
quast.py MF1_RagTagCorrScaff_TGSGapCloserONT.scaff_seqs -o TGS_Quast
```
# Preparing fasta for NCBI upload

## SeqKit
Before submitting a genome to NCBI terminal Ns need to be removed. I used code found on biostars: https://www.biostars.org/p/412636/
This was done using the below code:
```
seqkit -is replace -p "n+$" -r "" MF1_RagTagCorrScaff_TGSGapCloserONT.scaff_seqs
```

## Removing contigs smaller than 200 bp and rename contigs
original instructions
```
seqkit seq -m 200 MelFel_1.0.fa > MelFel_1.0_longcontigs.fa
for f in MelFel_1.0_longcontigs.fa; do
    seqkit replace -p '.+' -r 'contig{nr}' $f > $f.rename.fa
done
```
july 2024
```
seqkit seq -m 200 MelFal_2.0.fasta > MelFal_2.0_longcontigs.fasta
for f in MelFal_2.0_longcontigs.fasta; do
    seqkit replace -p '.+' -r 'contig{nr}' $f > $f.rename.fasta
done
```

## Convert multiline fasta to single line format
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < MelFel_1.0_longcontigs.fa.rename.fa > MelFel_1.0_longcontigs.fa.rename.oneline.fa
```
## Mapping to check for large contigs
```
samtools faidx MelFal_1.0.fa
bwa index MelFal_1.0.fa
bwa mem MelFal_1.0.fa MF1_S12_R1_001_val_1.fq.gz MF1_S12_R2_001_val_2.fq.gz > MF1.sam
samtools view -bh MF1.sam > MF1.bam
samtools sort -o MF1_sorted.bam MF1.bam
samtools index MF1_sorted.bam
qualimap bamqc --java-mem-size=8g -bam MF1_sorted.bam -outdir ./qualimap/nematode_qualimap -outfile MF1.pdf
```
Check which contigs are large and remove them by opening the file and deleting them. I removed Contigs: 180, 214, 226, 238, 240, 241, 243 from the fasta file as they had a coverage of over 500x

## Repeat modeler and masker
Really useful link: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.html#gsc.tab=0
### Repeat modeler
1. Build Database
```
conda activate repeatmodeler
BuildDatabase -name MelFal_postNCBI MelFal_1.0.fa
```
2. Run Repeatmodeler
```
module load miniconda
conda activate repeatmodeler
RepeatModeler -database MelFal_postNCBI -threads 28 -LTRStruct > run.out
```
### Repeat masker
```
RepeatMasker -lib MelFal_postNCBI-families.fa MelFal_1.0.fa
```
## Convert from a multiline fasta to a singleline fasta
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < MelFal_1.0.fa.masked > MelFal_1.0.masked.oneline.fa
```

## BWA mapping and Qualimap to calculate new coverage 
```
samtools faidx MelFal_1.0.masked.oneline.fa
bwa index MelFal_1.0.masked.oneline.fa
bwa mem MelFal_1.0.masked.oneline.fa MF1_S12_R1_001_val_1.fq.gz MF1_S12_R2_001_val_2.fq.gz > MF1.sam
samtools view -bh MF1.sam > MF1.bam
samtools sort -o MF1_sorted.bam MF1.bam
samtools index MF1_sorted.bam
qualimap bamqc --java-mem-size=8g -bam MF1_sorted.bam -outdir ./qualimap/nematode_qualimap -outfile MF1.pdf
```

## QUAST and BUSCO
```
quast.py <>.fasta -o FINAL_GENOME_Quast
busco -i <>.fasta -l nematoda_odb10 -o Mfallax_Meloidogyne -m genome --cpu 30
busco -i <>.fasta --auto-lineage -o Mfallax_AUTO -m genome --cpu 30
```

## FCS-GX 
NCBI use FCS-GX to check for contaminants (https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart)
This database is in the bioblob server location:
```
cd /nfs/bioblob/FCS
```
### Screen the genome
1. Assign the path to the --gx-db folder to GXDB_LOC.
```
GXDB_LOC=/nfs/bioblob/FCS/fcs_database
```
2. Retrieve the organism tax-id from NCBI Taxonomy: Fallax is 71801
3. Screen the genome:
```
python3 fcs.py screen genome --fasta MelFal_2.0_FINAL_1.fa --out-dir ./gx_out/ --gx-db "$GXDB_LOC" --tax-id 71801 
```
### Clean the genome
1. Perform cleaning actions on input genome:
```
zcat MelFal_2.0_FINAL_1.fa | python3 fcs.py clean genome --action-report ./gx_out/MelFal_2.0_FINAL_1.71801.fcs_gx_report.txt --output clean.fasta --contam-fasta-out contam.fasta
```
2. Split on internal contaminants instead of masking:
```
sed -i 's/FIX/SPLIT/g' ./gx_out/MelFal_2.0_FINAL_1.71801.fcs_gx_report.txt
zcat MelFal_2.0_FINAL_1.fa | python3 ./fcs.py clean genome --action-report ./gx_out/MelFal_2.0_FINAL_1.71801.fcs_gx_report.txt --output clean.fasta --contam-fasta-out contam.fasta
```

## Masking portions of the genome using Bedtools maskfasta
When I uploaded the genome to NCBI they said I had to mask the genome using Ns in certain areas. I used maskfasta in bedtools (https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html) I did this by:
```
bedtools maskfasta -fi MelFal_1.0_contigs_removed.fa -bed contig_removal.bed -fo MelFal_1.0_contigs_removed_bedtools.fa
```
