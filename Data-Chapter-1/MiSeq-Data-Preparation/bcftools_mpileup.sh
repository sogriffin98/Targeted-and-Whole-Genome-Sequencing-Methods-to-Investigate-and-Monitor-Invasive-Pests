#!/bin/bash

# Set variables
REFERENCE="/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/GCF_912470025.1_iVesVel2.1_genomic.fna"
BAM_DIR="/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/bam"
VCF_OUTPUT="/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/vcf/AH_48_indivs.vcf"

# Check if reference genome exists
if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference genome $REFERENCE not found!"
    exit 1
fi

# Check if BAM directory exists and contains BAM files
if [[ ! -d "$BAM_DIR" ]] || [[ -z $(ls "$BAM_DIR"/*primaryalignment.bam 2>/dev/null) ]]; then
    echo "Error: No primary alignment BAM files found in $BAM_DIR!"
    exit 1
fi

# Ensure output directory exists
mkdir -p $(dirname "$VCF_OUTPUT")

# Run bcftools command
bcftools mpileup -f "$REFERENCE" --annotate AD,DP --max-depth 10000 "$BAM_DIR"/*primaryalignment.bam \
    | bcftools call -m -v -f GQ --skip-variants indels -o "$VCF_OUTPUT"

# Check if the VCF file was created successfully
if [[ -f "$VCF_OUTPUT" ]]; then
    echo "VCF file successfully created: $VCF_OUTPUT"
else
    echo "Error: VCF file was not created!"
    exit 1
fi
