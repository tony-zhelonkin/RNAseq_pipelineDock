#!/bin/bash

# Define paths
BAM_DIR="/data/star_run/sep_run"  # Replace with the actual path to your BAM files
GTF_FILE="/data/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"  # Path to the annotation file

# Define the BAM file you want to skip
SKIP_FILE="AWe-01_S1_Aligned.sortedByCoord.out.bam"  # Replace with the actual filename you want to skip

# Loop through all BAM files in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Get the base name of the BAM file (without the directory and extension)
    BASE_NAME=$(basename "$BAM_FILE" _Aligned.sortedByCoord.out.bam)
    
    # Skip the specified file
    if [[ "$(basename "$BAM_FILE")" == "$SKIP_FILE" ]]; then
        echo "Skipping $BAM_FILE"
        continue
    fi
    
    # Run quant3p for the current BAM file
    quant3p -n "$BASE_NAME" -g "$GTF_FILE" "$BAM_FILE"
    
    # Optional: Check if quant3p ran successfully and log the result
    if [ $? -eq 0 ]; then
        echo "quant3p successfully processed $BAM_FILE"
    else
        echo "quant3p failed for $BAM_FILE" >&2
    fi
done

