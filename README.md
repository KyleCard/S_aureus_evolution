## Evolution under vancomycin selection drives divergent collateral sensitivity patterns in *Staphylococcus aureus*

Data, figures, and analysis code for the following manuscript:
[https://www.biorxiv.org/content/10.1101/2023.11.30.569373v3](https://www.biorxiv.org/content/10.1101/2023.11.30.569373v3)


## Overview

This repository contains:

- An R Notebook that loads curated variant calls and MIC measurements, performs statistical analyses, and renders the figures.
- Bash/Python pipeline scripts used to process short-read sequencing data with *breseq* and summarize read-level support for variants.
- MIC and curated variant files as input data for the R Notebook.
- Curated variant summaries and clustering results from the Bayesian latent-class analysis.
- Final figures as PDF/TIF.

If you only want to browse results, open the rendered notebook markdown: `S_aureus_evolution.md`.


## Repository layout

- `S_aureus_evolution.Rmd` — Primary analysis notebook (R Markdown). Produces the summary statistics and figures.
- `S_aureus_evolution.md` — Rendered markdown of the notebook for quick viewing.
- `data/` — Input and output data for/from the R Notebook.
	- `vancomycin_MICs.csv` — Vancomycin MIC time series for 18 evolving populations. **Input**
	- `collateral_drug_MICs.csv` — MICs for 8 antibiotics measured on ancestor and evolved lines (for collateral response analyses). **Input**
	- `mutations.csv` — Curated mutation calls (gene-level and metadata) used in the genomic analyses. **Input**
	- `variant_analysis/` — Outputs from read-level validation of variant calls:
		- `all_supporting_reads.tsv` — One row per read supporting a variant.
		- `variant_summary.tsv` — One row per variant per sample with coverage and support counts.
		- `variant_summary_curated.tsv` — Manually curated variant summary used downstream.
	- `cluster_results/` — Clustering assignments (K = 2 – 7) for control and vancomycin groups from BLCA consensus clustering. **Output**
- `figures/` — Final figures used in the manuscript.
	- Figures: `figure_1.pdf`, `figure_2.pdf`, `figure_3.pdf`, `figure_4.pdf`, `figure_S1.tif`, `figure_S2.tif`.
	- Cluster plots in `figures/cluster_plots/`.
- `sequencing_pipeline/` — Scripts to process raw reads and summarize variant support:
	- `1_trimmomatic.sh` — Quality trim paired-end reads.
	- `2_generate_updated_reference.sh` — Call ancestor vs ATCC_29213 and apply differences to create an updated reference.
	- `3_map_sequences.sh` — Run breseq (polymorphism mode) on experimental and control lines.
	- `4_identify_reads_supporting_variants.sh` — Extract BAM reads overlapping variant sites with samtools.
	- `5_summarize_base_qualities.py` — Parse supporting reads and compute per-variant support/quality summaries with pysam.


## Data dictionaries

1) `data/vancomycin_MICs.csv`
- `population` — Evolved population ID (1–18).
- `day` — Day sampled.
- `MIC` — Vancomycin MIC (µg/mL).

2) `data/collateral_drug_MICs.csv`
- `population` — Evolved population ID or `Ancestor` for the ancestral clone.
- `replicate` — Technical/biological replicate index.
- `antibiotic` — Drug abbreviation (see below).
- `MIC` — Minimum inhibitory concentration (µg/mL).
- `paired_ID` — Identifier linking ancestor/evolved measurements for response calculations.

Antibiotic abbreviations used in the notebook:
- `VAN` (vancomycin), `DAP` (daptomycin), `CFZ` (cefazolin), `SXT` (trimethoprim–sulfamethoxazole), `CLI` (clindamycin), `MEM` (meropenem), `NAF` (nafcillin), `GEN` (gentamicin).

3) `data/mutations.csv` (subset of key fields)
- `id` — Row identifier.
- `population` — Sample label (e.g., `T1`, `C37`).
- `DP` — Total read depth at the variant site.
- `AD` — Alternate-allele supporting read count.
- `median_phred` — Median base quality for variant-supporting bases.
- `median_error` — Approximate error rate from base qualities.
- `consensus_score` — *breseq* consensus evidence score.
- `polymorphism_score` — *breseq* polymorphism evidence score.
- `notes` — Manual curation notes; NA if kept; used to flag homopolymers, repeats, multicopy elements, and variants with low position-hash scores.
- `qualifying_mutation` — Boolean used to select variants for analysis.
- `vancomycin` — Indicator for treated or control group (per row).
- `evidence` — *breseq* evidence type (`RA` read-alignment, junctions for SVs, etc.).
- `seq_id`, `position` — Genomic coordinates (reference-dependent).
- `mutation_type` — SNP, small indel, intergenic, junction, etc.
- `mutation` — Human-readable change (e.g., `G -> T`, `(T)6 -> 5`, `H261N`).
- `freq` — Variant allele frequency (0–1 for polymorphisms, 1 for fixed).
- `annotation` — Functional annotation context (coding/intergenic and offsets).
- `gene` — Gene symbol used as a feature unit.
- `description` — Gene description.

4) `data/variant_analysis/all_supporting_reads.tsv`
- `sample` — Sample ID (e.g., `S1`, `C7`).
- `chrom`, `pos` — Genomic coordinates (1-based POS as in VCF).
- `ref`, `alt` — Reference and first ALT allele.
- `read` — Read name.
- `allele` — Observed allele in the read (`C`, `+INSSEQ`, `-LEN`).
- `phred` — Base quality at the variant position (if defined).

5) `data/variant_analysis/variant_summary.tsv`
- `sample`, `chrom`, `pos`, `ref`, `alt` — Variant key.
- `coverage` (DP) — Reads overlapping site.
- `alt_support` (AD) — Reads supporting first ALT allele.
- `median_phred` — Median base quality among supporting reads (if defined).

6) `data/variant_analysis/variant_summary_curated.tsv`
- Same as `variant_summary.tsv` with added `id` and `notes`, after manual curation.

7) `data/cluster_results/*.csv`
- `ID` — Sample (`T*` or `C*`).
- `Cluster` — Cluster assignment (e.g., `Cluster-2`) for the specified K.


## R Notebook summary (analysis outline)

The R Notebook (`S_aureus_evolution.Rmd`) performs the following:

1) Vancomycin resistance trajectories
- Loads `vancomycin_MICs.csv` and plots MIC over time for 18 evolving populations (Figure 1).

2) Genomic evolution and parallelism
- Loads `mutations.csv`, filters curated variants, constructs a population-by-gene mutation matrix, and computes Dice’s similarity coefficient within groups (treated vs control) with permutation testing (10,000 label shuffles) to assess significance (Figure S1).

3) Feature selection with weighted elastic-net logistic regression
- Uses *breseq* evidence scores to derive per-feature penalty factors (higher penalty for low-confidence variants) and fits elastic-net models with cross-validation (AUC). Bootstraps (1,000 iterations) to assess feature stability and derives odds ratios and CIs for selected genes (Figure 2B, 2C heatmap of selected features).

4) Bayesian latent class analysis (BLCA)
- Runs BLCA on selected gene sets within groups to identify mutation pattern classes and builds consensus matrices; performs hierarchical clustering and exports cluster assignments and heatmaps for K = 2–7 (figures saved to `figures/cluster_plots/`).

5) Collateral responses to first-line antibiotics
- Loads `collateral_drug_MICs.csv`, computes collateral response (log$_2$ relative MIC vs ancestor), summarizes by population and antibiotic, visualizes heatmaps/boxplots, and performs Mann–Whitney *U* tests and multiple-testing correction (Figure 3B–D).

6) Genotype–phenotype links and collateral response score (CSS)
- Tests associations between mutation pathways (e.g., WalKR regulon, *rpsU*, *yycH*) and collateral responses (point-biserial correlations, adjusted *p*-values), and estimates bootstrap CSS distributions by antibiotic and by cluster (ridgeline visualizations).

All figures are written to the `figures/` folder.


## How to run the notebook

Prerequisites (R packages, as used in the Rmd):

- CRAN: tidyverse, cowplot, proxy, foreach, doParallel, circlize, caret, glmnet, ggnewscale, BayesLCA, Matrix, RColorBrewer, ggsci, boot, ggridges
- Bioconductor: ComplexHeatmap

Install packages in R:

```r
install.packages(c(
	"tidyverse","cowplot","proxy","foreach","doParallel","circlize", "caret",
	"glmnet","ggnewscale","BayesLCA","Matrix","RColorBrewer","ggsci","boot","ggridges"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

Render the Notebook to GitHub-flavored markdown (Windows PowerShell example):

```powershell
# From the repository root
Rscript -e "rmarkdown::render('S_aureus_evolution.Rmd', output_format='rmarkdown::github_document')"
```

Notes:
- Parallel sections in the Rmd use all but one core by default; adjust if needed.
- Paths in the Rmd assume execution from the repository root.


## Sequencing pipeline (optional, for raw data reprocessing)

These scripts are designed for a Unix-like environment (Linux or Windows 11 with WSL2). They expect a directory structure with raw FASTQs under `raw_reads/` (ancestral, experimental_lines, control_lines). Outputs are written under `trimmed_read_files/` and `breseq_output/`.

Dependencies and versions used in scripts:
- Trimmomatic 0.39
- breseq 0.39.0 (includes `gdtools`)
- samtools 1.21
- Python 3.8+ with `pysam`

Example reprocessing steps:

1) Quality trim reads
```bash
bash sequencing_pipeline/1_trimmomatic.sh
```

2) Create an updated reference from the ancestor vs ATCC_29213
```bash
bash sequencing_pipeline/2_generate_updated_reference.sh
```

3) Map experimental and control lines with breseq (polymorphism mode)
```bash
bash sequencing_pipeline/3_map_sequences.sh
```

4) Extract reads overlapping variant sites
```bash
bash sequencing_pipeline/4_identify_reads_supporting_variants.sh
```

5) Summarize read-level support and qualities (produces TSVs in `data/variant_analysis/`)
```bash
python3 sequencing_pipeline/5_summarize_base_qualities.py
```

Manual curation: The repository includes a curated variant summary (`data/variant_analysis/variant_summary_curated.tsv`) derived from the automated summary. Downstream analyses use `data/mutations.csv`, which reflects this curation and additional annotations.


## Reproducing figures

- Run the R Markdown as above; figures are saved to `figures/` and subfolders.
- Pre-rendered versions are included:
	- Combined: `figures/figure_1.pdf`, `figure_2.pdf`, `figure_3.pdf`, `figure_4.pdf`.
	- Supplementary: `figures/figure_S1.tif`, `figure_S2.tif`.
	- Cluster heatmaps: `figures/cluster_plots/`.


## Quick start

If you want to:

- Browse results: open `S_aureus_evolution.md`.
- Inspect input data: see the `data/` folder and dictionaries above.
- Re-render analyses: install R packages and render `S_aureus_evolution.Rmd`.
- Reprocess raw reads: use the `sequencing_pipeline/` scripts on Linux/WSL.


## Citation

If you use this code or data, please cite the associated manuscript:

K. Card, D. Crozier, *et al.* (2025). Evolution under vancomycin selection drives divergent collateral sensitivity patterns in Staphylococcus aureus. [Manuscript].

An updated citation with DOI will be added upon publication. You may also cite this repository’s specific commit if referencing code.


## License and reuse

- All source code, notebooks, and documentation in this repository are licensed under the MIT License. See the LICENSE file at the repository root.
- Redistributed copies must retain the copyright notice and permission notice from the LICENSE file.
- Data and results in `data/` are provided under the same MIT terms.
- Third‑party tools, references, and embedded assets remain under their respective licenses.
- By submitting a contribution, you agree to license your contribution under the MIT License.


## Contact and support

- For questions, open a GitHub issue in this repository.
- Maintainers: Kyle Card (owner). Co-authors and contributors are acknowledged in the manuscript and notebook header.


