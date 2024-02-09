#!/bin/bash

# Demultiplexing, quality control, and adapter trimming of sequencing reads were performed with bcl-convert (v3.9.3) by SeqCenter.
# Here, we further quality filter the ancestral clone (Anc_R1_001 / Anc_R2_001) and evolved population (S*_R1_001 / S*_R2_001) 
# reads with Trimmomatic v0.39 and settings: SLIDINGWINDOW:4:20 MINLEN:36

# All paths are relative to this directory structure:

# .
# └── project/
#     ├── raw_reads
#     └── sequencing_pipeline/
#         ├── 1_trimmomatic.sh
#         ├── 2_generate_updated_reference.sh
#         └── 3_map_sequences.sh

# Create directories for trimmed read files
mkdir ../trimmed_read_files
mkdir ../trimmed_read_files/paired_forward_reads ../trimmed_read_files/paired_reverse_reads \
    ../trimmed_read_files/unpaired_forward_reads ../trimmed_read_files/unpaired_reverse_reads

# Trim reads
trimmomatic PE -threads 24 ../raw_reads/Anc_R1_001.fastq.gz ../raw_reads/Anc_R2_001.fastq.gz \
    ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz ../trimmed_read_files/unpaired_forward_reads/Anc_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz ../trimmed_read_files/unpaired_reverse_reads/Anc_R2_001.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:36

for i in {1..18}
do
    trimmomatic PE -threads 24 \
        ../raw_reads/S"$i"_R1_001.fastq.gz ../raw_reads/S"$i"_R2_001.fastq.gz \
        ../trimmed_read_files/paired_forward_reads/S"$i"_R1_001.fastq.gz ../trimmed_read_files/unpaired_forward_reads/S"$i"_R1_001.fastq.gz \
        ../trimmed_read_files/paired_reverse_reads/S"$i"_R2_001.fastq.gz ../trimmed_read_files/unpaired_reverse_reads/S"$i"_R2_001.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:36
done