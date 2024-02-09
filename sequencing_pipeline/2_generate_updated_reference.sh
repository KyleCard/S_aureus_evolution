#!/bin/bash

# Map the ancestral clone reads to the reference genome ATCC_29213.gbk with breseq 0.38.1 and settings: -j 24

# This step accounts for possible mutational differences between the ancestral clone used to found the replicate populations and the reference genome ATCC_29213 from NCBI. 
# These mutational differences will need to be discarded from the analysis of the evolved populations.

breseq -r ../ATCC_29213.gbk ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz -j 24 -o ../breseq_output/Ancestor

# Apply the mutational differences between the ancestral clone and the ATCC reference genome to the reference genome

gdtools APPLY -r ../ATCC_29213.gbk ../breseq_output/Ancestor/output/output.gd -f GFF3 -o ../Anc_updated.gff3

# Re-run breseq to verify

breseq -r ../Anc_updated.gff3 ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz -j 24 -o ../breseq_output/updated_ancestor