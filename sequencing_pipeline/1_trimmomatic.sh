#!/bin/bash

# Author: Kyle Card

# Here, we quality filter the ancestral clone (Anc_R1_001 / Anc_R2_001), evolved
# population (S*_R1_001 / S*_R2_001), and control lines (C*_R1_001 / C*_R2_001)
# paired-end Illumina sequencing reads from the raw_reads directory using Trimmomatic
# (v0.39) with the following parameters:
# - PE (paired-end reads)
# - threads: 24
# - SLIDINGWINDOW:4:20 (sliding window of 4 bases, cutting when average quality is below 20)
# - MINLEN:36 (discard reads shorter than 36 bases)

# All paths and code are relative to this directory structure:

#   .
#   └── project/
#       ├── raw_reads/
#       │   ├── control_lines/
#       │   │   ├── C1_R1_001.fastq.gz
#       │   │   ├── C1_R2_001.fastq.gz
#       │   │   └── etc.
#       │   ├── experimental_lines/
#       │   │   ├── S1_R1_001.fastq.gz
#       │   │   ├── S1_R2_001.fastq.gz
#       │   │   └── etc.
#       │   ├── Anc_R1_001.fastq.gz
#       │   └── Anc_R2_001.fastq.gz
#       └── sequencing_pipeline/
#           ├── 1_trimmomatic.sh
#           ├── 2_generate_updated_reference.sh
#           └── 3_map_sequences

## Create directory for trimmed read files
mkdir ../trimmed_read_files

# Create subdirectories for paired and unpaired reads
mkdir ../trimmed_read_files/paired_forward_reads \
    ../trimmed_read_files/paired_reverse_reads \
    ../trimmed_read_files/unpaired_forward_reads \
    ../trimmed_read_files/unpaired_reverse_reads

# Create subdirectories for control and experimental lines
mkdir ../trimmed_read_files/paired_forward_reads/control_lines \
    ../trimmed_read_files/paired_reverse_reads/control_lines \
    ../trimmed_read_files/unpaired_forward_reads/control_lines \
    ../trimmed_read_files/unpaired_reverse_reads/control_lines

mkdir ../trimmed_read_files/paired_forward_reads/experimental_lines \
    ../trimmed_read_files/paired_reverse_reads/experimental_lines \
    ../trimmed_read_files/unpaired_forward_reads/experimental_lines \
    ../trimmed_read_files/unpaired_reverse_reads/experimental_lines

# Trim ancestral reads
trimmomatic PE \
    -threads 24 \
    ../raw_reads/Anc_R1_001.fastq.gz \
    ../raw_reads/Anc_R2_001.fastq.gz \
    ../trimmed_read_files/paired_forward_reads/Anc_R1_001.fastq.gz \
    ../trimmed_read_files/unpaired_forward_reads/Anc_R1_001.fastq.gz \
    ../trimmed_read_files/paired_reverse_reads/Anc_R2_001.fastq.gz \
    ../trimmed_read_files/unpaired_reverse_reads/Anc_R2_001.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:36

# Trim experimental line reads
for i in {1..18}
do
    trimmomatic PE \
        -threads 24 \
        ../raw_reads/experimental_lines/S"$i"_R1_001.fastq.gz \
        ../raw_reads/experimental_lines/S"$i"_R2_001.fastq.gz \
        ../trimmed_read_files/paired_forward_reads/experimental_lines/S"$i"_R1_001.fastq.gz \
        ../trimmed_read_files/unpaired_forward_reads/experimental_lines/S"$i"_R1_001.fastq.gz \
        ../trimmed_read_files/paired_reverse_reads/experimental_lines/S"$i"_R2_001.fastq.gz \
        ../trimmed_read_files/unpaired_reverse_reads/experimental_lines/S"$i"_R2_001.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:36
done

# Trim control line reads
for i in {1..87}
do
    trimmomatic PE \
        -threads 24 \
        ../raw_reads/control_lines/C"$i"_R1_001.fastq.gz \
        ../raw_reads/control_lines/C"$i"_R2_001.fastq.gz \
        ../trimmed_read_files/paired_forward_reads/control_lines/C"$i"_R1_001.fastq.gz \
        ../trimmed_read_files/unpaired_forward_reads/control_lines/C"$i"_R1_001.fastq.gz \
        ../trimmed_read_files/paired_reverse_reads/control_lines/C"$i"_R2_001.fastq.gz \
        ../trimmed_read_files/unpaired_reverse_reads/control_lines/C"$i"_R2_001.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:36
done