#!/bin/bash

# Variables
GTF_FILE="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"  # Path to your GTF annotation file
INTERVAL_LIST_OUTPUT="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_ribosomal.interval_list"  # Output file for Picard

# Reference Genome Header (e.g., if you're using mouse genome GRCm39, update accordingly)
REFERENCE_GENOME="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"

# Extract rRNA from the GTF file and format it into Picard's interval_list format
echo "Generating ribosomal.interval_list from GTF file..."

# Add the header to the interval_list file
echo "@HD    VN:1.5  SO:coordinate" > "$INTERVAL_LIST_OUTPUT"

# Extract rRNA genes from the GTF and format them into interval_list format
grep -i "rRNA" "$GTF_FILE" | awk 'BEGIN {OFS="\t"} 
{
    if ($3 == "exon") {
        chr = $1;
        start = $4;
        end = $5;
        strand = $7;
        gene = $12;  # Update this based on your GTF file structure
        gsub(/"/, "", gene);
        gsub(/;/, "", gene);
        print chr, start, end, strand, gene;
    }
}' >> "$INTERVAL_LIST_OUTPUT"

echo "ribosomal.interval_list created: $INTERVAL_LIST_OUTPUT"

