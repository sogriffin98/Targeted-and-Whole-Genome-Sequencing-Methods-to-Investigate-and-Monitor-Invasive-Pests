#!/bin/bash

# Define input and output file patterns
INPUT_PATTERN="*_sorted.bam"
OUTPUT_SUFFIX=".mq30.maxinsert.primaryalignment.bam"

# Loop over all matching BAM files
for bam_file in $INPUT_PATTERN; do
    # Extract the sample name by removing the "_sorted.bam" suffix
    sample_name=$(basename "$bam_file" "_sorted.bam")

    # Define output filenames
    output_file="${sample_name}${OUTPUT_SUFFIX}"

    # Run bamtools filter
    bamtools filter -mapQuality '>=30' -isPrimaryAlignment 'true' -insertSize '<=800' -in "$bam_file" -out "$output_file"

    # Check if filtering was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed: $bam_file -> $output_file"

        # Index the filtered BAM file
        samtools index "$output_file"

        # Check if indexing was successful
        if [ $? -eq 0 ]; then
            echo "Successfully indexed: $output_file"
        else
            echo "Error indexing: $output_file"
        fi
    else
        echo "Error processing: $bam_file"
    fi
done
