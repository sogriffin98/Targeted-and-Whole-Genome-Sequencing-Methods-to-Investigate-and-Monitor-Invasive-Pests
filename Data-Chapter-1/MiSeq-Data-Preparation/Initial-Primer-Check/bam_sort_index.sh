#!/bin/bash

# Directory containing BAM files
BAM_DIR="./bam"

# Output directory (optional)
OUTPUT_DIR="./bam/sorted_bam"

# Loop through all BAM files in the directory
for bam_file in "$BAM_DIR"/*.bam; do
    if [ -f "$bam_file" ]; then
        # Extract the base filename without extension
        base_name=$(basename "$bam_file" .bam)
        
        echo "Processing $bam_file ..."
        
        # Sort the BAM file
        sorted_bam="$OUTPUT_DIR/${base_name}_sorted.bam"
        samtools sort -o "$sorted_bam" "$bam_file"
        
        # Index the sorted BAM file
        samtools index "$sorted_bam"
        
        echo "Sorted and indexed: $sorted_bam"
    fi
done

echo "All BAM files processed."
