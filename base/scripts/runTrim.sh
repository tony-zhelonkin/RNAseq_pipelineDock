#!/bin/bash

# Number of parallel FastQC jobs
JOB_NUMBER=6  # Number of parallel FastQC processes to run

# Define the patterns for R1 and R2
R1_PATTERNS=("R1" "_1" "_R1_001" "_R1_paired")
R2_PATTERNS=("R2" "_2" "_R2_001" "_R2_paired")

# Trimming parameters
LEADING=3
TRAILING=3
SLIDINGWINDOW="4:20"
MINLEN=36
ILLUMINACLIP_SETTINGS="2:30:10"
ADAPTERS="/usr/share/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# Directories for trimmed reads and QC results
TRIMMED_DIR="./trimmed/reads"
FASTQC_OUTPUT_DIR="./trimmed/postTrimQC/fastQC"
MULTIQC_OUTPUT_DIR="./trimmed/postTrimQC/multiQC"
LOG_DIR="./trimmed/logs"

# Create output directories if they don't exist
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${FASTQC_OUTPUT_DIR}"
mkdir -p "${MULTIQC_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Loop through all FASTQ files and process paired-end reads
for R1_PATTERN in "${R1_PATTERNS[@]}"; do
  for R2_PATTERN in "${R2_PATTERNS[@]}"; do
    # Find files that match the R1 pattern
    for R1_FILE in *${R1_PATTERN}*.fastq*; do
      # Generate the corresponding R2 filename
      R2_FILE=$(echo "${R1_FILE}" | sed "s/${R1_PATTERN}/${R2_PATTERN}/")

      # Check if the R2 file exists
      if [[ -f "${R2_FILE}" ]]; then
        # Define output file names without unpaired reads
        SAMPLE_NAME=$(basename "${R1_FILE}" | sed "s/${R1_PATTERN}.*//")
        TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_paired_trimmed.fastq.gz"
        TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_paired_trimmed.fastq.gz"
        LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_trimmomatic.log"

        # Run Trimmomatic (no unpaired output files)
        echo "Processing ${SAMPLE_NAME} - R1: ${R1_FILE}, R2: ${R2_FILE}"
        trimmomatic PE -threads 8 \
          "${R1_FILE}" "${R2_FILE}" \
          "${TRIMMED_R1}" /dev/null \
          "${TRIMMED_R2}" /dev/null \
          ILLUMINACLIP:"${ADAPTERS}":"${ILLUMINACLIP_SETTINGS}" \
          LEADING:"${LEADING}" TRAILING:"${TRAILING}" \
          SLIDINGWINDOW:"${SLIDINGWINDOW}" MINLEN:"${MINLEN}" \
          > "${LOG_FILE}" 2>&1

        # Add trimmed files to a list for FastQC processing
        echo "${TRIMMED_R1}" >> "${LOG_DIR}/fastqc_files_list.txt"
        echo "${TRIMMED_R2}" >> "${LOG_DIR}/fastqc_files_list.txt"
      fi
    done
  done
done

# Function to run FastQC and log output
run_fastqc() {
    local file=$1
    local base_name=$(basename "$file")
    local log_file="${LOG_DIR}/${base_name}_fastqc.log"

    echo "Running FastQC on $file..."
    fastqc "$file" --outdir "${FASTQC_OUTPUT_DIR}" --quiet > "$log_file" 2>&1
}

export -f run_fastqc
export FASTQC_OUTPUT_DIR LOG_DIR

# Run FastQC on each trimmed file in parallel using GNU parallel with logging
cat "${LOG_DIR}/fastqc_files_list.txt" | parallel -j "$JOB_NUMBER" run_fastqc {}

# Run MultiQC to aggregate the FastQC reports
echo "Running MultiQC on FastQC results"
multiqc "${FASTQC_OUTPUT_DIR}" -o "${MULTIQC_OUTPUT_DIR}"

echo "Trimming and QC complete."
