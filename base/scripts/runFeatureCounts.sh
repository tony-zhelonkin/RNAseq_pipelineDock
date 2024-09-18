#!/bin/bash

# Directory containing BAM files
BAM_DIR="."  # Adjust this if your BAM files are in a different directory

# Annotation file
ANNOTATION_FILE="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"

# Output directories
RAW_OUTPUT_DIR="./featureCounts/raw"
MATRIX_OUTPUT_DIR="./featureCounts/count_matrices"
LOG_DIR="./featureCounts/logs"

mkdir -p "$RAW_OUTPUT_DIR"
mkdir -p "$MATRIX_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Number of threads for featureCounts
THREADS=12  # Adjust this based on your system's capabilities

# Output files
COUNTS_FILE="$RAW_OUTPUT_DIR/counts_matrix.txt"
SORTED_COUNTS_FILE="$MATRIX_OUTPUT_DIR/sorted_counts_matrix.txt"
LOG_FILE="$LOG_DIR/featureCounts.log"

# Find all BAM files in the directory and sort them lexicographically
mapfile -t SORTED_BAM_FILES < <(find "$BAM_DIR" -type f -name "*.bam" | sort)

echo "Processing BAM files: ${SORTED_BAM_FILES[*]}"

# Run featureCounts on all sorted BAM files in one command
featureCounts \
    -a "$ANNOTATION_FILE" \
    -o "$COUNTS_FILE" \
    -p \
    -s 0 \
    -M \
    --fraction \
    -T "$THREADS" \
    "${SORTED_BAM_FILES[@]}" \
    > "$LOG_FILE" 2>&1

# Reorder columns in the output file
awk 'BEGIN {FS=OFS="\t"} 
    NR==1 {next}  # Skip the first line (featureCounts summary)
    NR==2 {
        printf "Geneid";
        for (i=7; i<=NF; i++) {
            split($i, a, "/");  # Split the path
            split(a[length(a)], b, ".");  # Split the filename
            printf "\t%s", b[1];  # Print only the sample name
        }
        printf "\n";
        next;
    }
    {
        printf "%s", $1;
        for (i=7; i<=NF; i++) printf "\t%s", $i;
        printf "\n";
    }' "$COUNTS_FILE" > "$SORTED_COUNTS_FILE"

# Optional: Logging completion
echo "featureCounts analysis and count matrix creation completed." >> "$LOG_FILE"
echo "Completed processing all BAM files. Results saved to $SORTED_COUNTS_FILE."
