#!/usr/bin/env python3

# Author: Kyle Card

"""
In the previous step, we used samtools to identify reads from BAM files
that support breseq variant calls. In this step, we use the pysam library
to parse the alignment for each supporting read, pinpointing the exact position
within that read that corresponds to the variant base call. Once we have that
position, we will extract the corresponding base quality score from the quality
string in the BAM record for that read. These steps will allow us to
calculate the median base quality score for the variant base calls across all
experimental and control lines.
"""

"""
Extract all reads that support the ALT allele in every breseq-processed
sample and write

1  all_supporting_reads.tsv   (one row per supporting read)
2  variant_summary.tsv        (one row per variant, incl. zero-support cases)

Directory layout assumed:

breseq_output/
  experimental_lines/S1/data/…
  experimental_lines/S2/data/…
  …
  control_lines/C1/data/…
  …

Each “data” folder must contain:
    reference.bam   – coordinate-sorted, indexed
    output.vcf      – breseq VCF
    identified_reads/supporting_reads.sorted.bam
                     – result of the earlier step, indexed
"""

import csv, sys
from pathlib import Path
from statistics import median
import pysam

ROOT   = Path("../breseq_output")          # edit if needed
GROUPS = [("experimental_lines", 18, "S"),
          ("control_lines",      87, "C")]


def load_vcf(vcf_path):
    """Return {(chrom,pos): dict} with REF/ALT and variant type."""
    variants = {}
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, *_ = line.rstrip().split("\t")
            pos  = int(pos)
            alts = alt.split(",")
            entry = {"ref": ref, "alts": alts}

            if all(len(a) == len(ref) == 1 for a in alts):
                entry["type"] = "snp"
            elif all(len(a) > len(ref) for a in alts):
                entry["type"]    = "ins"
                entry["alt_ins"] = [a[1:] for a in alts]   # drop anchor
            elif all(len(a) < len(ref) for a in alts):
                entry["type"]    = "del"
                entry["del_len"] = len(ref) - min(len(a) for a in alts)
            else:
                entry["type"] = "complex"

            variants[(chrom, pos)] = entry
    return variants


def check_read(rec, pos, var):
    """Return (allele, phred) if read supports the FIRST ALT allele; else None."""
    alt1 = var["alts"][0]

    if var["type"] == "snp":
        rp = rec.get_reference_positions()
        if pos - 1 in rp:
            q_idx = rp.index(pos - 1)
            base  = rec.query_sequence[q_idx]
            if base == alt1:
                return base, rec.query_qualities[q_idx]

    elif var["type"] == "ins":
        pairs  = rec.get_aligned_pairs()          # indels included by default
        anchor = [p for p in pairs if p[1] == pos - 1 and p[0] is not None]
        if not anchor:
            return
        q_anchor = anchor[-1][0]
        ins = []
        i = q_anchor + 1
        while i < len(rec.query_sequence) and pairs[i][1] is None:
            ins.append(rec.query_sequence[i]); i += 1
        if "".join(ins) == var["alt_ins"][0]:
            quals = rec.query_qualities[q_anchor + 1 : q_anchor + 1 + len(ins)]
            return f"+{''.join(ins)}", min(quals)        # min Phred over insertion

    elif var["type"] == "del":
        pairs = rec.get_aligned_pairs()
        gap   = [p for p in pairs
                 if pos - 1 <= (p[1] or pos) < pos - 1 + var["del_len"]
                 and p[0] is None]
        if len(gap) == var["del_len"]:
            return f"-{var['del_len']}", None            # Phred undefined for deletion

    return None


reads_f   = open("all_supporting_reads.tsv", "w", newline="")
summary_f = open("variant_summary.tsv",      "w", newline="")
reads_w   = csv.writer(reads_f,   delimiter="\t")
summary_w = csv.writer(summary_f, delimiter="\t")

reads_w.writerow(["sample","chrom","pos","ref","alt",
                  "read","allele","phred"])
summary_w.writerow(["sample","chrom","pos","ref","alt",
                    "coverage","alt_support","median_phred"])


for subdir, n, prefix in GROUPS:
    for i in range(1, n + 1):
        sample = f"{prefix}{i}"
        data   = ROOT / subdir / sample / "data"
        bam_fp = data / "identified_reads/supporting_reads.sorted.bam"
        vcf_fp = data / "output.vcf"

        if not bam_fp.exists():
            print(f"# WARNING: missing {bam_fp}", file=sys.stderr)
            continue

        variants = load_vcf(vcf_fp)
        bam      = pysam.AlignmentFile(bam_fp, "rb")

        for (chrom, pos), var in variants.items():
            alt1        = var["alts"][0]          # first ALT only
            coverage    = bam.count(chrom, pos - 1, pos)      # all reads
            alt_support = 0
            phreds      = []

            if var["type"] != "complex":
                for rec in bam.fetch(chrom, pos - 1, pos):
                    res = check_read(rec, pos, var)
                    if res:
                        allele, phred = res
                        reads_w.writerow([sample, chrom, pos,
                                          var["ref"], alt1,
                                          rec.query_name, allele, phred])
                        alt_support += 1
                        if phred is not None:
                            phreds.append(phred)

            median_phred = "." if not phreds else median(phreds)
            summary_w.writerow([sample, chrom, pos,
                                var["ref"], alt1,
                                coverage, alt_support, median_phred])

        bam.close()

reads_f.close()
summary_f.close()

print("Done.\n  • all_supporting_reads.tsv\n  • variant_summary.tsv")


"""
We then manually curated the variant summary file, generating the
variant_summary_curated.tsv file.

"""