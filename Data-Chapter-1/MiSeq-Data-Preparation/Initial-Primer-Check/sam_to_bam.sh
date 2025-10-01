#!/bin/bash

OUTPUT_DIR="/nfs/sgriffin_data2/MiSeq_Run_1/250124_M08033_0096_000000000-GP2VP/Alignment_1/20250125_083656/Fastq/bwa/bam"

# Loop through all .sam files in the current directory
for sam_file in *.sam; do
    # Ensure the file exists before processing
    if [[ -f "$sam_file" ]]; then
        # Define the output BAM file name
        bam_file="${OUTPUT_DIR}/$(basename "$sam_file" .sam).bam"
        
        echo "Converting $sam_file to $bam_file..."
        
        # Convert SAM to BAM
        samtools view -Sb "$sam_file" > "$bam_file"
        
        if [[ $? -eq 0 ]]; then
            echo "Successfully converted $sam_file to $bam_file"
        else
            echo "Error converting $sam_file"
        fi
    fi
done

echo "Conversion completed!"

