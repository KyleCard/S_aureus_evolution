#!/bin/bash

# Map the evolved population reads to the updated reference genome (Anc_updated.gff3) with breseq 0.38.1 and settings: -j 24 -p

for i in {1..18}
do
    breseq -r ../Anc_updated.gff3 ../trimmed_read_files/paired_forward_reads/S"$i"_R1_001.fastq.gz ../trimmed_read_files/paired_reverse_reads/S"$i"_R2_001.fastq.gz -j 24 -p -o ../breseq_output/S"$i"
done