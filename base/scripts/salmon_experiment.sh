#!/bin/bash

# Number of parallel jobs (for GNU parallel)
JOB_NUMBER=6  # Number of parallel Salmon processes to run

# Directories (use absolute paths)
FASTQ_DIR="."  # Directory where FASTQ files are located
SALMON_OUTPUT_DIR="./salmon"  # Output directory for Salmon results
LOG_DIR="./salmon_logs"  # Directory for logs

# Create output directories if they don't exist
mkdir -p "$SALMON_OUTPUT_DIR" "$LOG_DIR"

# Prompt user for Salmon index information
read -p "Do you want to build a new Salmon index? (yes/no): " build_index

if [ "$build_index" == "yes" ]; then
    # Prompt user for the path to the FASTA file for indexing
    read -p "Enter the path to the FASTA file for building the index: " fasta_file

    # Prompt user for the name of the index directory
    read -p "Enter the name of the directory to store the Salmon index: " index_name

    # Create the index directory
    index_dir="$SALMON_OUTPUT_DIR/$index_name"

    # Build the Salmon index
    echo "Building Salmon index in $index_dir..."
    salmon index -t "$fasta_file" -i "$index_dir"

    # Update SALMON_INDEX to point to the newly built index
    SALMON_INDEX="$index_dir"
elif [ "$build_index" == "no" ]; then
    # Prompt user for the path to the existing Salmon index
    read -p "Enter the path to the existing Salmon index directory: " SALMON_INDEX
else
    echo "Invalid input. Exiting."
    exit 1
fi

# Function to infer strandedness using Salmon
run_salmon() {
    local fastq_file=$1
    local base_name=$(basename "$fastq_file")
    local sample_name="${base_name%.*}"  # Remove file extension for sample name
    local result_file="$SALMON_OUTPUT_DIR/${sample_name}_salmon_output.txt"
    local log_file="$LOG_DIR/${sample_name}_salmon.log"

    echo "Running Salmon on $fastq_file..."
    salmon quant -i "$SALMON_INDEX" -l A -r "$fastq_file" -o "$result_file" > "$log_file" 2>&1
}

export -f run_salmon
export SALMON_OUTPUT_DIR LOG_DIR SALMON_INDEX

# Find all FASTQ files and run Salmon in parallel with logging
find "$FASTQ_DIR" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fasta" \) | parallel -j "$JOB_NUMBER" run_salmon {}

# Optional: Logging completion
echo "Salmon analysis completed." >> "$LOG_DIR/salmon_analysis.log"
