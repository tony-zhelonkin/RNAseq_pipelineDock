#!/bin/bash

# Number of parallel jobs (for GNU parallel)
JOB_NUMBER=6  # Number of parallel Picard, RSeQC, and Samtools processes to run

# Directories (use absolute paths)
BAM_DIR="."  # Directory containing BAM files
PICARD_OUTPUT_DIR="../post_align_QC/picard"  # Output directory for Picard metrics
RSEQ_OUTPUT_DIR="../post_align_QC/rseqc"  # Output directory for RSeQC metrics
SAMTOOLS_OUTPUT_DIR="../post_align_QC/samtools"  # Output directory for Samtools flagstat
MULTIQC_DIR="../post_align_QC/multiQC"  # Output directory for MultiQC report
LOG_DIR="../post_align_QC/logs"  # Directory to store logs

# Define the reference genome and refFlat file (make sure these paths are correct)
REFERENCE_GENOME="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"  # Path to reference genome
REF_FLAT_FILE="/mouse_genome/ref_genome/GRCm39_refFlat"  # Path to refFlat file for Picard
BED_ANNOTATION_FILE="/mouse_genome/ref_genome/GRCm39.bed"  # BED file for RSeQC infer_experiment.py

# Create output directories if they don't exist
mkdir -p "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# Function to run Picard CollectRnaSeqMetrics
run_picard() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local log_file="$LOG_DIR/${base_name}_picard.log"
    echo "Running Picard CollectRnaSeqMetrics on $bam_file..."
    picard CollectRnaSeqMetrics \
        I="$bam_file" \
        O="$PICARD_OUTPUT_DIR/${base_name}_rna_metrics.txt" \
        R="$REFERENCE_GENOME" \
        REF_FLAT="$REF_FLAT_FILE" \
        STRAND_SPECIFICITY="FIRST_READ_TRANSCRIPTION_STRAND" \
        RIBOSOMAL_INTERVALS="ribosomal.interval_list" \
        > "$log_file" 2>&1
}

# Function to run RSeQC infer_experiment.py
run_rseqc() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local log_file="$LOG_DIR/${base_name}_rseqc.log"
    echo "Running RSeQC infer_experiment.py on $bam_file..."
    infer_experiment.py \
        -i "$bam_file" \
        -r "$BED_ANNOTATION_FILE" \
       # -q 30 \
        -s 200000 \
        > "$log_file" 2>&1
}

# Function to run Samtools flagstat
run_samtools_flagstat() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local output_file="$SAMTOOLS_OUTPUT_DIR/${base_name}_flagstat.txt"
    local log_file="$LOG_DIR/${base_name}_samtools_flagstat.log"
    echo "Running Samtools flagstat on $bam_file..."
    samtools flagstat "$bam_file" > "$output_file" 2> "$log_file"
}

# Export the functions and variables so they can be used by GNU parallel
export -f run_picard
export -f run_rseqc
export -f run_samtools_flagstat
export PICARD_OUTPUT_DIR RSEQ_OUTPUT_DIR SAMTOOLS_OUTPUT_DIR REFERENCE_GENOME REF_FLAT_FILE BED_ANNOTATION_FILE LOG_DIR

# Find all BAM files and run Picard, RSeQC, and Samtools flagstat in parallel
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_picard {}
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_rseqc {}
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_samtools_flagstat {}

# Run MultiQC to aggregate the results
echo "Running MultiQC..."
multiqc "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" -o "$MULTIQC_DIR" -n "PostAlignmentQC"

echo "Picard, RSeQC, Samtools flagstat, and MultiQC processing complete."
