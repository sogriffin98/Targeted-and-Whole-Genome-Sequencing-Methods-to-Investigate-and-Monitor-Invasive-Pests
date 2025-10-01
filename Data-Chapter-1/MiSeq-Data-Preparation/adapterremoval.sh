#!/bin/bash

# Set directories
input_dir="/nfs/sgriffin_data2/MiSeq_Run_2/250205_M08033_0097_000000000-GP4D8/Alignment_1/20250206_082834/Fastq"
output_dir="nfs/sgriffin_data2/MiSeq_Run_2/250205_M08033_0097_000000000-GP4D8/Alignment_1/20250206_082834/Fastq/AdapterRemoval"

# Loop through each sample (paired-end reads)
for file in "$input_dir"/*_R1_001.fastq; do
    # Check if the file exists
    if [ ! -e "$file" ]; then
        echo "No R1 files found in $input_dir"
        exit 1
    fi

    # Extract the sample name (assuming file naming convention includes '_R1_001.fastq')
    sample_name=$(basename "$file" "_R1_001.fastq")

    # Check if paired-end file exists
    file_R2="$input_dir/${sample_name}_R2_001.fastq"
    if [ ! -e "$file_R2" ]; then
        echo "No paired-end file found for $sample_name. Skipping..."
        continue
    fi

    # Define output prefix
    output_prefix="$output_dir/$sample_name"

    # Run AdapterRemoval for paired-end data
    echo "Processing paired-end sample: $sample_name"
    AdapterRemoval --file1 "$file" \
                   --file2 "$file_R2" \
                   --basename "$output_prefix" \
                   --trimns \
                   --trimqualities \
                   --collapse

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $sample_name"
    else
        echo "Error processing $sample_name"
    fi
done
