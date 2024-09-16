#!/bin/bash

# Start timing the script
START_TIME=$(date +%s)

# Define paths
GENOME_DIR="/mouse_genome/ref_index100"
BATCH_1_DIR="."
OUTPUT_DIR="./STAR_alignment/bam"
LOG_DIR="./STAR_alignment/logs"

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Define naming conventions for read pairs
R1_PATTERNS=("_R1" "_0001" "_1")
R2_PATTERNS=("_R2" "_0002" "_2")

# Function to log and output a command (only show errors in terminal)
log_and_run() {
    local CMD="$1"
    local LOG_FILE="$2"

    echo "Running command and logging to $LOG_FILE"  # Minimal terminal output
    eval "$CMD" >> "$LOG_FILE" 2>> >(tee -a "$LOG_FILE" >&2)  # Log everything, but show errors in terminal
}

# Function to process a pair of FASTQ files
process_fastq_pair() {
    local R1_FILE="$1"
    local R2_FILE="$2"
    local BASE_NAME="$3"
    local BATCH_DIR="$4"

    # Define read group (@RG) information
    local ID="${BASE_NAME}"
    local SM="Sample_${BASE_NAME}"
    local PL="ILLUMINA"
    local PU="${BATCH_DIR##*/}_${BASE_NAME}"

    # Create log file for this sample
    local LOG_FILE="${LOG_DIR}/${BASE_NAME}_log.txt"

    echo "Processing sample: $BASE_NAME"  # Minimal terminal output

    # STAR alignment
    log_and_run "STAR --runThreadN 4 \
                 --genomeDir \"$GENOME_DIR\" \
                 --readFilesIn \"$R1_FILE\" \"$R2_FILE\" \
                 --readFilesCommand zcat \
                 --sjdbOverhang 100 \
                 --limitBAMsortRAM 25000000000 \
                 --outFileNamePrefix \"$OUTPUT_DIR/${BASE_NAME}_\" \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
                 --twopassMode Basic" "$LOG_FILE"

    echo "Finished processing sample: $BASE_NAME"  # Minimal terminal output
}

# Export function and variables for parallel execution
export -f process_fastq_pair log_and_run
export GENOME_DIR OUTPUT_DIR LOG_DIR

# Identify R1 files based on naming patterns
is_r1_file() {
    local FILE="$1"
    for pattern in "${R1_PATTERNS[@]}"; do
        if [[ "$FILE" =~ $pattern ]]; then
            echo "true"
            return
        fi
    done
    echo "false"
}

# Find the corresponding R2 file
find_r2_file() {
    local R1_FILE="$1"
    for i in "${!R1_PATTERNS[@]}"; do
        local R2_FILE="${R1_FILE/${R1_PATTERNS[$i]}/${R2_PATTERNS[$i]}}"
        if [[ -f "$R2_FILE" ]]; then
            echo "$R2_FILE"
            return
        fi
    done
    return 1  # No matching R2 file found
}

# Main loop: Find all R1 files, process them in parallel
find "$BATCH_1_DIR" -type f -name "*.fastq.gz" -o -name "*.fq.gz" | while read -r FILE; do
    if [[ "$(is_r1_file "$FILE")" == "true" ]]; then
        BASE_NAME=$(basename "$FILE" | sed -E 's/(_R1|_0001|_1)(_paired)?(_001)?\.(fastq|fq)\.gz//')

        R2_FILE=$(find_r2_file "$FILE")
        if [[ -n "$R2_FILE" ]]; then
            BATCH_DIR=$(dirname "$FILE")
            echo "$FILE $R2_FILE $BASE_NAME $BATCH_DIR"
        else
            echo "Warning: Corresponding R2 file not found for $BASE_NAME" >&2
        fi
    fi
done | parallel --progress -j 2 --memfree 10G --colsep ' ' process_fastq_pair {1} {2} {3} {4}

# End timing and calculate total duration
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# Save total time to log file
echo "Total time: $TOTAL_TIME seconds" | tee "$LOG_DIR/time.txt"
