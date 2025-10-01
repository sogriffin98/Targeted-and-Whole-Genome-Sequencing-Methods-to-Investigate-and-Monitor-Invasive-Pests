#!/bin/bash

conda activate bwa

# Define reference genome
REFERENCE="/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/GCF_912470025.1_iVesVel2.1_genomic.fna"

# Number of threads
THREADS=4

# Index the reference genome
samtools faidx "$REFERENCE"
bwa index "$REFERENCE"

# Loop over paired-end FASTQ files
for READ in /nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/AdapterRemoval/*.fastq; do
    SAMPLE=$(basename "$READ" .fastq)
    
    # Align reads using bwa aln
    bwa aln -t "$THREADS" "$REFERENCE" "$READ" > "/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/sam/${SAMPLE}.sai"
    
    # Convert SAI to SAM using bwa samse
    bwa samse "$REFERENCE" "/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/sam/${SAMPLE}.sai" "$READ" > "/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/sam/${SAMPLE}.sam"
done
