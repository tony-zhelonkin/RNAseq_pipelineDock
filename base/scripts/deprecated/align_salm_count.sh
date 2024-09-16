#!/bin/bash

# Define paths
GENOME_DIR="/data/ref_re_index"  # Replace with actual path to STAR index
GTF_FILE="/data/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"  # Replace with actual GTF file path
BATCH_1_DIR="/data/Ada_Macrophages_bulkRNAseq_2024/Run_1"  # Replace with actual path to batch_1 FASTQ files
BATCH_2_DIR="/data/Ada_Macrophages_bulkRNAseq_2024/Run_2"  # Replace with actual path to batch_2 FASTQ files
OUTPUT_DIR="/data/star_run/sep_run"  # Replace with actual path for outputs

# Load modules or activate environment if necessary (e.g., conda activate your_env)

# Process each batch
for BATCH_DIR in "$BATCH_1_DIR" "$BATCH_2_DIR"; do
    # Loop through FASTQ files (assuming paired-end reads with _R1.fastq.gz and _R2.fastq.gz naming convention)
    for R1_FILE in "$BATCH_DIR"/*_R1_001.fastq.gz; do
        # Get base name and corresponding R2 file
        BASE_NAME=$(basename "$R1_FILE" _R1_001.fastq.gz)
        R2_FILE="$BATCH_DIR/${BASE_NAME}_R2_001.fastq.gz"
        
        # Define read group (@RG) information
        ID="${BASE_NAME}"  # Use the base name as ID
        SM="Sample_${BASE_NAME}"  # Sample name
        PL="ILLUMINA"  # Platform
        PU="${BATCH_DIR##*/}_${BASE_NAME}"  # Platform unit, includes batch and base name
        
        # STAR alignment
        STAR --runThreadN 8 \
             --genomeDir "$GENOME_DIR" \
             --readFilesIn "$R1_FILE" "$R2_FILE" \
             --readFilesCommand zcat \
	     --sjdbOverhang 107 \
	     --limitBAMsortRAM 36000000000 \
             --outFileNamePrefix "$OUTPUT_DIR/${BASE_NAME}_" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
             --twopassMode Basic

        # Infer strandedness with Salmon
        SALMON_OUTPUT_DIR="${OUTPUT_DIR}/${BASE_NAME}_salmon"
        salmon quant -i "$GENOME_DIR" \
                     -l A \
                     -1 "$R1_FILE" \
                     -2 "$R2_FILE" \
                     --gcBias \
                     --validateMappings \
                     -p 8 \
                     -o "$SALMON_OUTPUT_DIR"
                     
        # Determine inferred strandedness from Salmon logs
        STRANDEDNESS=$(grep "observed library type" "$SALMON_OUTPUT_DIR/logs/salmon_quant.log" | awk '{print $NF}')
        echo "Inferred strandedness for $BASE_NAME: $STRANDEDNESS"

        # Convert Salmon's inferred strandedness to HTSeq-count parameter
        if [[ "$STRANDEDNESS" == "ISR" || "$STRANDEDNESS" == "FR" ]]; then
            STRANDED_OPTION="yes"
        elif [[ "$STRANDEDNESS" == "ISF" || "$STRANDEDNESS" == "RF" ]]; then
            STRANDED_OPTION="reverse"
        else
            STRANDED_OPTION="no"
        fi

        # Run htseq-count on the resulting BAM file
        BAM_FILE="${OUTPUT_DIR}/${BASE_NAME}_Aligned.sortedByCoord.out.bam"
        htseq-count -f bam -r pos -s "$STRANDED_OPTION" -t exon -i gene_id "$BAM_FILE" "$GTF_FILE" > "${OUTPUT_DIR}/${BASE_NAME}_htseq_counts.txt"
    done
done

