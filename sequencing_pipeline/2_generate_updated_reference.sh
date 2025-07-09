#!/bin/bash

# Author: Kyle Card

# Map the ancestral clone reads to the reference genome ATCC_29213.gbk
# with breseq 0.39 and settings:
# -j: number of threads to use (24)
# -o: output directory for breseq results
# -r: reference genome in GenBank format

# This step accounts for possible mutational differences between the
# ancestral clone used to found the replicate populations and the
# reference genome ATCC_29213 from NCBI. These mutational differences
# will need to be discarded from the analysis of the evolved populations.

breseq \
    -j 24 \
    -o ../breseq_output/ancestor \
    -r ../ATCC_29213.gbk \
    ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz

# Apply the mutational differences between the ancestral clone and the ATCC
# reference genome to the reference genome using gdtools and the following
# settings:
# -r: reference genome in GenBank format
# -f: output format (GFF3)
# -o: output file for the updated reference genome in GFF3 format

gdtools \
    APPLY \
    -f GFF3 \
    -o ../Anc_updated.gff3 \
    -r ../ATCC_29213.gbk \
    ../breseq_output/ancestor/output/output.gd \

# Re-run breseq to verify

breseq \
    -j 24 \
    -o ../breseq_output/updated_ancestor \
    -r ../Anc_updated.gff3 \
    ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz