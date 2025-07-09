#!/bin/bash

# Author: Kyle Card

# This script uses samtools (v 1.21) to identify reads that support variants.

# For a given variant position specified in a VCF file, we first identify all the
# reads in the BAM file that cover this position and show the alternative allele.

# The BAM and VCF files are outputted by the breseq pipeline and are located in the
# ../breseq_output/experimental_lines/S"$i"/data directory for experimental lines and
# ../breseq_output/control_lines/C"$i"/data for control lines.

# The VCF file is named "output.vcf" and the BAM file is named "reference.bam".

set -euo pipefail

extract_reads() {
    local sample_dir=$1                      # …/data directory
    local bam="$sample_dir/reference.bam"
    local vcf="$sample_dir/output.vcf"
    local outdir="$sample_dir/identified_reads"
    mkdir -p "$outdir"

    # 1. Variant positions → BED (0‑based start, 1‑based end)
    awk -v OFS='\t' 'BEGIN{FS="\t"} !/^#/ {print $1, $2-1, $2}' "$vcf" > "$outdir/vars.bed"

    # 2. Grab every read that overlaps ≥1 variant position
    samtools view -h -L "$outdir/vars.bed" -b "$bam" > "$outdir/supporting_reads.bam"

    # 3. Coordinate‑sort (if needed) and index
    samtools sort -o "$outdir/supporting_reads.sorted.bam" "$outdir/supporting_reads.bam"
    samtools index "$outdir/supporting_reads.sorted.bam"
}

# Experimental lines S1–S18
for i in {1..18}; do
    extract_reads "../breseq_output/experimental_lines/S${i}/data"
done

# Control lines C1–C87
for i in {1..87}; do
    extract_reads "../breseq_output/control_lines/C${i}/data"
done