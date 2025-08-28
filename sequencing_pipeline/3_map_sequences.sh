#!/bin/bash

# Author: Kyle Card

# Map the evolved population reads to the updated reference genome
# (Anc_updated.gff3) with breseq 0.39 and settings:
# -j: number of threads to use (24)
# -p: polymorphism mode
# -o: output directory for breseq results
# -r: updated reference genome in GFF3 format

# Experimental lines (VAN_1 to VAN_18)
for i in {1..18}
do
    breseq \
    -j 24 \
    -p \
    -o ../breseq_output/experimental_lines/VAN_"$i" \
    -r ../Anc_updated.gff3 \
    ../trimmed_read_files/paired_forward_reads/experimental_lines/VAN_"$i"_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/experimental_lines/VAN_"$i"_R2_001.fastq.gz
done

# Control lines (CTL_1 to CTL_87)
for i in {1..87}
do
    breseq \
    -j 24 \
    -p \
    -o ../breseq_output/control_lines/CTL_"$i" \
    -r ../Anc_updated.gff3 \
    ../trimmed_read_files/paired_forward_reads/control_lines/CTL_"$i"_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/control_lines/CTL_"$i"_R2_001.fastq.gz
done