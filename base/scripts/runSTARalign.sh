#!/bin/bash

# Start timing the script
START_TIME=$(date +%s)

# Define paths
GENOME_DIR="/mouse_genome/ref_index100"
BATCH_1_DIR="."
OUTPUT_DIR="./STAR_alignment/bam"
LOG_DIR="./STAR_alignment/logs"

# Define computational resources
JOB_NUMBER="1" # number of GNU parallel jobs
GB_FREE="10G" # number of GB left free before any new process gets spawned

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

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
    log_and_run "STAR --runThreadN 10 \
                 --genomeDir \"$GENOME_DIR\" \
                 --readFilesIn \"$R1_FILE\" \"$R2_FILE\" \
                 --readFilesCommand zcat \
                 --sjdbOverhang 100 \
                 --limitBAMsortRAM 45000000000 \
                 --outFileNamePrefix \"$OUTPUT_DIR/${BASE_NAME}_\" \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
                 --twopassMode Basic" "$LOG_FILE"

    echo "Finished processing sample: $BASE_NAME"  # Minimal terminal output
}

# Export function and variables for parallel execution
export -f process_fastq_pair log_and_run
export GENOME_DIR OUTPUT_DIR LOG_DIR

# Identify if a file is an R1 file
is_r1_file() {
    local FILE="$1"
    if [[ "$FILE" =~ (_R1|_1)\.(fastq|fq|fasta)\.gz$ ]]; then
        echo "true"
    else
        echo "false"
    fi
}

# Generate the corresponding R2 file from an R1 file
get_r2_file() {
    local R1_FILE="$1"
    if [[ "$R1_FILE" =~ _R1 ]]; then
        echo "${R1_FILE/_R1/_R2}"
    elif [[ "$R1_FILE" =~ _1 ]]; then
        echo "${R1_FILE/_1/_2}"
    else
        return 1  # No valid R1 pattern found
    fi
}

# Main loop: Find all R1 files, process them in parallel
find "$BATCH_1_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fasta.gz" \) | while read -r FILE; do
    if [[ "$(is_r1_file "$FILE")" == "true" ]]; then
        BASE_NAME=$(basename "$FILE" | sed -E 's/(_R1|_1)(_paired)?(_001)?\.(fastq|fq|fasta)\.gz//')

        # Generate corresponding R2 file name
        R2_FILE=$(get_r2_file "$FILE")
        
        # Check if R2 file exists
        if [[ -f "$R2_FILE" ]]; then
            BATCH_DIR=$(dirname "$FILE")
            echo "$FILE $R2_FILE $BASE_NAME $BATCH_DIR"
        else
            echo "Warning: Corresponding R2 file not found for $BASE_NAME" >&2
        fi
    fi
done | parallel --progress -j "$JOB_NUMBER" --memfree "$GB_FREE" --colsep ' ' process_fastq_pair {1} {2} {3} {4}

# End timing and calculate total duration
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# Save total time to log file
echo "Total time: $TOTAL_TIME seconds" | tee "$LOG_DIR/time.txt"
