# Evolution under vancomycin selection drives divergent collateral sensitivity patterns in <i>Staphylococcus aureus</i>

**Authors**: Kyle Card and Arda Durmaz

**Last revised**: 8/9/2025

**Description**: In this study, we investigate how <i>Staphylococcus
aureus</i> evolution under vancomycin selection affects its
susceptibilities to several first-line antibiotics. We examined the
genetic basis of these changes using whole-genome sequencing and
identified key mutations associated with altered collateral responses.
Our findings provide insights into the mechanisms underlying collateral
sensitivity and highlight the potential for exploiting these patterns in
clinical settings.

The following code implements the analysis pipeline described in the
study.

**Citation**: K. Card, D. Crozier, *et al.* (2025). Evolution under
vancomycin selection drives divergent collateral sensitivity patterns in
<i>Staphylococcus aureus</i>. *bioRxiv*, doi:
<https://doi.org/10.1101/2023.11.30.569373>.

## Prerequisites - required packages, file paths, and parameters

``` r
library(tidyverse) # For data manipulation (dplyr, tidyr) and plotting (ggplot2)
library(cowplot) # For combining ggplot2 plots and themes
library(proxy) # For calculating pairwise distances and similarities
library(foreach) # For parallel processing
library(doParallel) # For parallel processing support
library(circlize) # For circular visualization and color scales

# Install ComplexHeatmap package with BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap) # For creating complex heatmaps

library(caret) # For creating cross-validation folds
library(glmnet) # For fitting elastic-net models
library(ggnewscale) # For adding multiple color scales to ggplot2 plots
library(BayesLCA) # For Bayesian latent class analysis
library(Matrix) # For sparse matrix operations
library(RColorBrewer) # For color palettes
library(ggsci) # For scientific color palettes
library(boot) # For bootstrapping functions
library(ggridges) # For ridge plots
```

``` r
# Suppress summarise info messages from dplyr
options(dplyr.summarise.inform = FALSE)
```

``` r
# File paths for input data and schema
input_vanc_MICs <- "data/vancomycin_MICs.csv"
input_experimental_schema <- "figures/PDFs/experimental_schema.pdf"
input_mutations <- "data/mutations.csv"
input_collateral_MICs <- "data/collateral_drug_MICs.csv"

# File paths for output data and plots
output_cv_plot <- "plots/cv_performance_plot.png"
output_feature_freq <- "data/feature_stability.tsv"
output_final_heatmap <- "plots/final_mutation_heatmap.pdf"
```

``` r
CONFIDENCE_FREQ_THRESHOLD <- 0.95 # VAF threshold to distinguish consensus from
# polymorphism
MIN_FEATURE_COUNT <- 3 # Minimum number of samples a feature must appear in
POSITIVE_CLASS_WEIGHT <- 5 # Weight for the positive class in the model
PERMUTATIONS <- 10000 # Number of permutations for similarity test
CV_FOLDS <- 5 # Number of folds for cross-validation
BOOTSTRAP_ITERATIONS <- 1000 # Number of bootstrap samples for feature stability
MIN_FEATURE_STABILITY <- 0.2 # Minimum stability frequency for features to be
# included in the final model
BLCA_CLUSTERS <- 5 # Number of clusters for BLCA
BLCA_ITERATIONS <- 1000 # Number of iterations for BLCA
CRS_ITERATIONS <- 100 # Number of iterations for CRS bootstrapping
```

------------------------------------------------------------------------

# Evolution of vancomycin-intermediate resistance in MSSA

The following data represent the minimum inhibitory concentrations
(MICs) of the 18 experimental populations evolved under increasing
vancomycin concentrations until they reached intermediate resistance
levels (4 - 8 $\mu$g/mL).

``` r
vancomycin_MICs <- read_csv(
  file = input_vanc_MICs,
  show_col_types = FALSE
)

print(vancomycin_MICs)
```

    ## # A tibble: 234 × 3
    ##    population   day MIC  
    ##         <dbl> <dbl> <chr>
    ##  1          1     1 1    
    ##  2          2     1 1.5  
    ##  3          3     1 1.5  
    ##  4          4     1 1    
    ##  5          5     1 1.5  
    ##  6          6     1 1.5  
    ##  7          7     1 1.5  
    ##  8          8     1 1    
    ##  9          9     1 1    
    ## 10         10     1 1.5  
    ## # ℹ 224 more rows

``` r
# --- Function that plots vancomycin MIC data with the experimental schema
# (Figure 1) ---

generate_figure_1B <- function(.x) {
  # In the MIC column, replace all "-" with "NA"
  .x$MIC <- .x$MIC %>%
    str_replace_all(
      pattern = "-",
      replacement = "NA"
    ) %>%
    as.numeric()

  # Remove all rows with NA values in the MIC column
  .x <- .x %>%
    filter(!is.na(MIC))

  # Plot of vancomycin MICs over time
  vanc_plot <- .x %>%
    ggplot(
      aes(
        x = day,
        y = MIC
      )
    ) +
    geom_rect(
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = 4,
        ymax = 8
      ),
      fill = "#E1EEF4",
      alpha = 0.1
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    facet_wrap(
      ~population,
      nrow = 3,
      ncol = 6
    ) +
    scale_x_continuous(breaks = seq(1, 14, 2)) +
    scale_y_continuous(
      limits = c(0, 8),
      breaks = seq(0, 8, 2)
    ) +
    xlab("Day") +
    ylab(
      expression("MIC (\u03BCg/mL)")
    ) +
    theme_cowplot() +
    theme(
      panel.grid.major = element_line(
        color = "#3a3a3a",
        linewidth = 0.2
      ),
      plot.background = element_rect(fill = "#FFFFFF")
    )

  return(vanc_plot)
}

figure_1B <- generate_figure_1B(.x = vancomycin_MICs)

print(figure_1B)
```

![](S_aureus_evolution_files/figure-gfm/Generate%20Figure%201-1.png)<!-- -->

------------------------------------------------------------------------

# Populations follow distinct adaptive pathways under vancomycin selection

We used the `breseq` pipeline (v0.39.0) to identify mutations in the
control and treated populations. The control lines evolved in TSB medium
without vancomycin, whereas the treated lines evolved in TSB medium
supplemented with vancomycin. We concatenated the `breseq` output files
for each control and treated population, labeled “C1” through “C87” and
“T1” through “T18”, respectively. We also manually curated the `breseq`
output files to ensure that the mutations were correctly identified and
annotated. There were a total of **4,936** variants.

We excluded **33** mutations within multicopy elements, such as
ribosomal RNA operons and synthetases responsible for charging tRNA.
Gene conversions may cause these mutations, but short-read sequencing
data cannot fully resolve them.

For structural variants, including indels less than 50 bp, large
deletions and duplications identified by new junction (JC) evidence, we
relied on a distinct set of quality metrics that assess the integrity of
split-read alignments. A putative junction was considered a true
positive only if it satisfied both of the following criteria:

1.  The junction had to be supported by reads that aligned to different
    starting positions in the reference sequence, as indicated by a
    position-hash score of 3 or higher. This criterion filters out
    artifacts arising from PCR.

2.  The junction had to be supported by reads mapping to both strands,
    thereby mitigating artifacts common to a single strand that occur
    during library preparation or sequencing. All putative structural
    variants passing these filters were visually inspected to confirm
    the supporting read alignments. If a variant did not satisfy either
    criterion, we excluded it as a false positive. **33** structural
    variants were excluded in this step.

Additionally, we manually curated SNPs and 1-bp indels supported by read
alignment (RA) evidence. If the mutation had a low position-hash score,
or occurred in a homopolymer tract or a tandem repeat region, we
classified it as a false positive. **169** mutations were excluded in
this step.

A total of **4,701** mutations passed the above quality control steps.

Moreover, following Deatherage et al. 2017 and Card et al. 2021, for
each remaining mutation, we determined whether it *qualified* for
further analyses based on the following criteria:

1.  It was a nonsynonymous point mutation or small indel in a single
    gene.
2.  An intergenic mutation within 150 bp upstream of the start of a
    gene.
3.  A large deletion if at least one of the affected genes was also
    found to be mutated in another population.

Conversely, we excluded mutations from further analyses if they were:

1.  Synonymous (does not affect the resulting amino acid sequence).
2.  Multicopy elements (e.g., ribosomal RNA operons) that may result
    from gene conversions but cannot be fully resolved using short-read
    sequencing data.

A total of **3,703** mutations qualify based on these criteria: **200**
mutations in the vancomycin-evolved lines and **3,503** mutations in the
control lines.

``` r
# --- Load and preprocess data ---

variants <- read_csv(
  file = input_mutations,
  show_col_types = FALSE
)

print(variants)
```

    ## # A tibble: 4,942 × 20
    ##       id population    DP    AD median_phred median_error consensus_score
    ##    <dbl> <chr>      <dbl> <dbl>        <dbl>        <dbl>           <dbl>
    ##  1     1 T1           142   141           34       0.0004            748.
    ##  2     2 T1           122   122           34       0.0004            511 
    ##  3     3 T1           109    10           34       0.0004            373.
    ##  4     4 T1           157     8           34       0.0004            619.
    ##  5     5 T1           133   133           34       0.0004            581.
    ##  6     6 T1           158     8           34       0.0004            617.
    ##  7     7 T1           167   167           34       0.0004            734.
    ##  8     8 T1           191    10           34       0.0004            758.
    ##  9     9 T1           164    12           34       0.0004            614.
    ## 10    10 T1           176    12           34       0.0004            634.
    ## # ℹ 4,932 more rows
    ## # ℹ 13 more variables: polymorphism_score <dbl>, notes <chr>,
    ## #   qualifying_mutation <lgl>, vancomycin <lgl>, evidence <chr>, seq_id <dbl>,
    ## #   position <dbl>, mutation_type <chr>, mutation <chr>, freq <dbl>,
    ## #   annotation <chr>, gene <chr>, description <chr>

``` r
variants_processed <- variants %>%
  # Filter out likely false positives based on manual data inspection. We made
  # notes on whether the variant occurred in a tandem repeat, homopolymer tract,
  # or multicopy element. All other variants were given NA.
  filter(is.na(notes)) %>%
  # Ensure score columns are numeric and handle non-finite values
  mutate(
    across(
      .cols = c(consensus_score, polymorphism_score),
      ~ {
        numeric_values <- as.numeric(.)
        if_else(
          condition = numeric_values < 0,
          true = NA_real_,
          false = numeric_values
        )
      }
    ),
    # Replace infinite polymorphism scores with the max finite score
    polymorphism_score = if_else(
      is.infinite(polymorphism_score),
      max(polymorphism_score[is.finite(polymorphism_score)], na.rm = TRUE),
      polymorphism_score
    )
  ) %>%
  # Calculate scaled weights for regularization. Low-confidence variants get a
  # higher weight, meaning they will be penalized more.
  mutate(
    score_type = case_when(
      evidence != "RA" ~ "indel",
      as.numeric(freq) <= CONFIDENCE_FREQ_THRESHOLD ~ "polymorphism",
      as.numeric(freq) > CONFIDENCE_FREQ_THRESHOLD ~ "consensus",
      TRUE ~ "other"
    ),
    score_to_use = case_when(
      score_type == "consensus" ~ consensus_score,
      score_type == "polymorphism" ~ polymorphism_score,
      TRUE ~ NA_real_
    )
  ) %>%
  # Group by score type to calculate median and median absolute deviation (MAD)
  # for each distribution separately
  mutate(
    .by = score_type,
    local_m = median(
      x = score_to_use,
      na.rm = TRUE
    ),
    local_mad = mad(
      x = score_to_use,
      na.rm = TRUE
    ),
    score_scaled = case_when(
      score_type == "indel" ~ 1.0,
      is.na(score_to_use) ~ 10.0, # High penalty for missing scores
      score_to_use <= (local_m - 2 * local_mad) ~ 10.0, # High penalty for low scores
      score_to_use <= (local_m + 2 * local_mad) ~ 1.0, # Medium penalty
      score_to_use > (local_m + 2 * local_mad) ~ 0.1, # Low penalty for high scores
      TRUE ~ 1.0 # Default case
    )
  ) %>%
  # Filter to qualifying mutations and create a clean variant ID
  filter(qualifying_mutation == TRUE) %>%
  rename(variant_id = gene)

# --- Create Mutation and Weight Matrices ---

# Create the binary mutation matrix representing the qualifying gene-level
# mutations in each evolutionary replicate
mutation_matrix_long <- variants_processed %>%
  distinct(population, variant_id) %>%
  mutate(value = 1)

mutation_matrix <- mutation_matrix_long %>%
  pivot_wider(
    names_from = variant_id,
    values_from = value,
    values_fill = 0
  )

# Create the weight matrix (average scaled score for each variant in each
# population)
weight_matrix_long <- variants_processed %>%
  group_by(population, variant_id) %>%
  summarise(
    mean_score = mean(
      x = score_scaled,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

weight_matrix <- weight_matrix_long %>%
  pivot_wider(
    names_from = variant_id,
    values_from = mean_score
  )


# --- Align and Finalize Matrices ---

# Ensure rows and columns are in the same order for both matrices.
sample_order <- mutation_matrix$population
feature_order <- setdiff(
  x = colnames(mutation_matrix),
  y = "population"
)

# Convert to numeric matrices and remove the population column.
# We use this matrix for the genomic similarity analysis.
mutation_matrix <- as.matrix(mutation_matrix[, feature_order])
rownames(mutation_matrix) <- sample_order

W_temp <- as.matrix(weight_matrix[, feature_order])
rownames(W_temp) <- weight_matrix$population

# Align W_temp to X and fill missing weights with a default value of 1.
W <- W_temp[sample_order, feature_order]
W[is.na(W)] <- 1

# Filter features that don't appear in a minimum number of samples
feature_counts <- colSums(mutation_matrix)
features_to_keep <- which(feature_counts >= MIN_FEATURE_COUNT)

mutation_matrix_final <- mutation_matrix[, features_to_keep]
W_final <- W[, features_to_keep]
final_feature_names <- colnames(mutation_matrix_final)

# Define the response variable (y) and observation weights
response_var <- as.factor(
  ifelse(
    str_starts(
      sample_order, "T"
    ),
    yes = "Treated",
    no = "Control"
  )
)

observation_weights <- ifelse(
  response_var == "Treated",
  yes = POSITIVE_CLASS_WEIGHT,
  no = 1.0
)
```

We quantified the extent of genomic parallelism in the
vancomycin-adapted lines and in the control lines by calculating Dice’s
Similarity Coefficient (*S*) for each population pair, where:

$$S = \frac{2|X \cap Y|}{|X|+|Y|}$$

$|X|$ and $|Y|$ represent the number of genes with qualifying mutations
in each population, and $|X \cap Y|$ is the number of mutated genes in
common between them. *S* therefore ranges from 0, when the pair share no
mutations in common, to 1, when both have mutations in exactly the same
set of genes (Deatherage et al., 2017; Card et al. 2021).

The perform_genomic_analysis function will output a list containing two
elements:

1.  Average pairwise similarity within the vancomycin-treated group
    $S_v \approx 0.063$ and within the control group
    ($S_c \approx 0.023$). In other words, two populations that
    independently evolved under vancomycin selection had, on average,
    **6.3%** of their mutated genes in common, whereas those that
    evolved under identical conditions, but without vancomycin, shared
    on average only **2.3%** of their mutated genes.

2.  The difference between these averages ($S_v - S_c \approx 0.04$)
    represents how much greater (on average) the vancomycin-adapted
    lines’ similarity is relative to the control lines. We expect that
    the vancomycin-adapted lines will exhibit a higher degree of genomic
    parallelism than the control lines.

3.  To evaluate the statistical significance of the observed difference
    between the vancomycin-adapted and control lines, the code generates
    a distribution under the null hypothesis that there is no difference
    between the two groups. To achieve this distribution, we permuted
    the population labels by resampling without replacement 10,000
    times. For each permutation, we calculated the difference in $S_v$
    and $S_c$. We then compared our observed difference to the null
    distribution of differences. The *p*-value is the proportion of
    permutations where the difference was greater than or equal to the
    observed difference.

``` r
# --- Function to calculate pairwise similarity and average similarity within
# groups. If permutation = TRUE, shuffles the population labels before
# calculating similarity. ---

perform_genomic_analysis <- function(.x, num_perm = 10000) {
  calculate_similarity <- function(.x, permutation = FALSE) {
    # Sample the lines without replacement (shuffles population labels)
    if (permutation == TRUE) {
      shuffled_row_names <- sample(row.names(.x))
      row.names(.x) <- shuffled_row_names
    }

    # Compute Dice's Similarity Coefficient for all possible replicate pairs.
    pwise_simil_matrix <- as.matrix(
      simil(
        x = .x,
        method = "Dice",
        by_rows = TRUE,
        upper = TRUE,
        diag = TRUE
      )
    )

    # All values above (and including) the matrix diagonal are converted to NA
    diag(pwise_simil_matrix) <- NA
    pwise_simil_matrix[upper.tri(pwise_simil_matrix)] <- NA

    # Convert matrix to a long data frame for analysis
    pwise_simil_df <- pwise_simil_matrix %>%
      as.data.frame() %>%
      rownames_to_column(var = "population_1") %>%
      pivot_longer(
        cols = -population_1,
        names_to = "population_2",
        values_to = "value",
        values_drop_na = TRUE
      ) %>%
      mutate(
        group_1 = str_remove_all(string = population_1, pattern = "\\d"),
        group_2 = str_remove_all(string = population_2, pattern = "\\d"),
        comparison_type = ifelse(
          group_1 == group_2,
          yes = "within_group",
          no = "between_group"
        )
      )

    # Computes the mean pairwise similarity within each group
    avg_similarity_df <- pwise_simil_df %>%
      filter(comparison_type == "within_group") %>%
      group_by(group_1) %>%
      summarize(
        avg = mean(value),
        .groups = "drop"
      )

    # Calculate the difference between groups
    control_avg <- avg_similarity_df %>%
      filter(group_1 == "C") %>%
      pull(avg)

    treated_avg <- avg_similarity_df %>%
      filter(group_1 == "T") %>%
      pull(avg)

    difference <- treated_avg - control_avg
    names(difference) <- "difference"

    if (permutation == FALSE) {
      results <- list(avg_similarity_df, difference)
    } else {
      results <- avg_similarity_df %>%
        mutate(group = ifelse(
          group_1 == "C",
          yes = "control",
          no = "treated"
        )) %>%
        dplyr::select(-group_1) %>%
        pivot_wider(
          names_from = group,
          values_from = avg
        ) %>%
        mutate(difference = treated - control)
    }
    return(results)
  }

  # --- Function to perform permutation test in parallel ---

  perform_permutation_test <- function(.x, num_perm) {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    results <- foreach(
      i = seq_len(num_perm),
      .combine = rbind,
      .packages = c("tidyverse", "proxy"),
      .export = "calculate_similarity"
    ) %dopar% {
      set.seed(i)
      calculate_similarity(
        .x = .x,
        permutation = TRUE
      )
    }
    stopCluster(cl)
    return(results)
  }

  # --- Result outputs ---

  # Calculate observed similarity
  similarity <- calculate_similarity(.x = .x)

  # Perform permutation test
  permutation_results <- perform_permutation_test(
    .x = .x,
    num_perm = num_perm
  )

  # Calculate significance
  observed_difference <- similarity[[2]]
  significance <- sum(
    permutation_results$difference >= observed_difference
  ) / num_perm

  return(
    list(
      similarity = similarity,
      permutation_results = permutation_results,
      significance = significance
    )
  )
}


# --- Run the genomic analysis ---

genomic_analysis <- perform_genomic_analysis(
  .x = mutation_matrix,
  num_perm = PERMUTATIONS
)


# --- Results ---

similarity <- genomic_analysis$similarity
names(similarity) <- c("Average pairwise similarity", "Difference in means")

genomic_results <- list(
  similarity_results = similarity,
  permutation_results = genomic_analysis$permutation_results,
  p_value = genomic_analysis$significance
)

print(genomic_results)
```

    ## $similarity_results
    ## $similarity_results$`Average pairwise similarity`
    ## # A tibble: 2 × 2
    ##   group_1    avg
    ##   <chr>    <dbl>
    ## 1 C       0.0233
    ## 2 T       0.0631
    ## 
    ## $similarity_results$`Difference in means`
    ## difference 
    ## 0.03981136 
    ## 
    ## 
    ## $permutation_results
    ## # A tibble: 10,000 × 3
    ##    control treated difference
    ##      <dbl>   <dbl>      <dbl>
    ##  1  0.0198  0.0238   0.00394 
    ##  2  0.0208  0.0166  -0.00420 
    ##  3  0.0214  0.0137  -0.00764 
    ##  4  0.0206  0.0175  -0.00309 
    ##  5  0.0207  0.0204  -0.000291
    ##  6  0.0203  0.0200  -0.000281
    ##  7  0.0205  0.0227   0.00213 
    ##  8  0.0207  0.0217   0.000989
    ##  9  0.0212  0.0162  -0.00497 
    ## 10  0.0203  0.0186  -0.00168 
    ## # ℹ 9,990 more rows
    ## 
    ## $p_value
    ## [1] 0

``` r
# --- Plot the distribution of differences in mean similarity (Figure S1) ---

figure_S1 <- genomic_analysis$permutation_results %>%
  ggplot(
    aes(x = difference)
  ) +
  geom_histogram(
    binwidth = 0.001,
    fill = "#7C8DA2"
  ) +
  geom_vline(
    xintercept = similarity[[2]],
    linetype = "dashed",
    color = "#000000"
  ) +
  scale_x_continuous(
    limits = c(-0.02, 0.04),
    breaks = seq(-0.02, 0.04, by = 0.01)
  ) +
  annotate(
    geom = "text",
    x = similarity[[2]] - 0.01,
    y = 850,
    label = "Observed difference",
    color = "#000000"
  ) +
  annotate(
    geom = "text",
    x = similarity[[2]] - 0.012,
    y = 800,
    label = "P",
    color = "#000000",
    fontface = "italic"
  ) +
  annotate(
    geom = "text",
    x = similarity[[2]] - 0.008,
    y = 800,
    label = "< 0.0001",
    color = "#000000"
  ) +
  xlab("Difference in mean similarity") +
  ylab("Frequency") +
  theme_cowplot() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white")
  )

print(figure_S1)
```

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](S_aureus_evolution_files/figure-gfm/Plot%20the%20permuted%20distribution%20of%20differences-1.png)<!-- -->

The permutation test showed that the vancomycin-adapted lines differ
significantly in genomic parallelism relative to the control lines.
However, this analysis does not identify which genes are associated with
vancomycin selection. To identify these genes, we implemented a weighted
elastic-net logistic regression using the `glmnet` package. This
approach incorporates both L1 (LASSO) and L2 (Ridge) regularization to
perform feature selection on a high-dimensional dataset. We included
qualifying mutations that occurred in at least three populations,
resulting in **484** gene-level features for the model.

A key feature of our model is the integration of variant-calling
confidence scores to guide the regularization process. Each gene-level
feature was assigned a specific penalty weight based on its underlying
`breseq` confidence scores. This process was implemented in a
discretized manner, where variants with low, medium, and high confidence
were assigned regularization weights of 10.0, 1.0, and 0.1,
respectively. This weighting scheme forces the model to apply a stronger
penalty to low-confidence variants, thereby prioritizing the selection
of features supported by high-quality data.

To identify a robust set of informative features, we first determined
the optimal model hyperparameters ($\alpha$ and $\lambda$) via k-fold
cross-validation. We then performed a bootstrap analysis over 1,000
iterations, fitting the weighted model with these optimal parameters to
a new bootstrap sample of the data in each run. Given the
sparsity-inducing property of L1 regularization, we tracked the
frequency of non-zero coefficients for each gene across all iterations.
Mutations exhibiting non-zero coefficients in at least 20% of the
iterations were retained, constituting the final set of stable features
for subsequent Bayesian latent class analysis (BLCA). Finally, we
computed the 95% confidence intervals for each selected gene’s
coefficient and exponentiated the endpoints to obtain the odds ratio,
allowing us to determine the direction and strength of its association
with the selection condition.

``` r
# Visualization of how different variants will be penalized.
# Red indicates a higher weight and thus a stronger penalty.
weight_col_fun <- colorRamp2(
  breaks = c(0, 5, 10),
  colors = c("white", "orange", "firebrick")
)

# Create heatmap of weights for each population and feature
weight_heatmap <- Heatmap(
  W_final,
  name = "Scaled Penalty",
  col = weight_col_fun,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 4),
  border_gp = gpar(col = "grey"),
  row_split = response_var, # Split columns (samples) by group
  show_row_dend = FALSE,
  show_column_dend = FALSE, # Also hide column dendrogram for clarity
  row_title = "Populations",
  column_names_rot = 45
)

draw(weight_heatmap)
```

![](S_aureus_evolution_files/figure-gfm/Heatmap%20of%20regularization%20weights-1.png)<!-- -->

``` r
# --- Cross-validation to find optimal hyperparameters ---

# Create stratified folds to ensure class balance in each fold
set.seed(123) # for reproducibility

cv_folds_list <- createFolds(
  response_var,
  k = CV_FOLDS,
  list = TRUE,
  returnTrain = FALSE
)

# Convert the list of folds into a vector of fold IDs for each observation
cv_fold_ids <- numeric(length(response_var))

# Assign fold numbers to each observation based on the created folds
for (i in 1:length(cv_folds_list)) {
  cv_fold_ids[cv_folds_list[[i]]] <- i
}

# Define the grid of alpha values to test
alpha_grid <- c(0.25, 0.5, 0.75, 1.0)

# Run CV for each alpha
cv_results <- map_df(alpha_grid, function(current_alpha) {
  # cv.glmnet performs k-fold cross-validation to find the optimal lambda
  cv_fit <- cv.glmnet(
    x = mutation_matrix_final,
    y = response_var,
    family = "binomial",
    alpha = current_alpha,
    penalty.factor = colMeans(W_final), # Penalty factors are per-feature
    weights = observation_weights,
    foldid = cv_fold_ids, # Use the correctly formatted fold assignments
    type.measure = "auc" # Use AUC for performance evaluation
  )

  # Store results
  tibble(
    alpha = current_alpha,
    lambda_min = cv_fit$lambda.min, # Lambda giving max AUC
    lambda_1se = cv_fit$lambda.1se, # Simplest model within 1-SE of max AUC
    auc_min = max(cv_fit$cvm),
    n_coef_min = cv_fit$nzero[which.max(cv_fit$cvm)],
    auc_1se = cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.1se],
    n_coef_1se = cv_fit$nzero[cv_fit$lambda == cv_fit$lambda.1se]
  )
})

# Choose the best parameters
best_params <- cv_results %>%
  arrange(desc(auc_min)) %>%
  head(1)
chosen_alpha <- best_params$alpha
chosen_lambda <- best_params$lambda_min


# --- Bootstrapping ---

# To assess how robustly features are selected, we train the model on many
# bootstrap samples of the data and record how often each feature's
# coefficient is non-zero.

# List to store coefficients from each bootstrap iteration
bootstrap_coefs <- map(1:BOOTSTRAP_ITERATIONS, function(i) {
  # Create a bootstrap sample (sampling with replacement)
  set.seed(123 + i)
  sample_indices <- sample(1:nrow(mutation_matrix_final), replace = TRUE)

  # Fit the model with the chosen parameters on the bootstrap sample
  fit <- glmnet(
    x = mutation_matrix_final[sample_indices, ],
    y = response_var[sample_indices],
    family = "binomial",
    alpha = chosen_alpha,
    lambda = chosen_lambda,
    penalty.factor = colMeans(W_final),
    weights = observation_weights[sample_indices]
  )

  # Return the full vector of coefficients
  as.vector(coef(fit))
})

# Convert the list of coefficient vectors into a matrix
coef_matrix <- do.call(rbind, bootstrap_coefs)

# Remove the first column corresponding to the intercept
coef_matrix <- coef_matrix[, -1]

colnames(coef_matrix) <- colnames(mutation_matrix_final)


# --- Calculate stability and odds ratios from bootstrap distributions ---

# Iterate over each feature (column) and calculate stats
odds_ratio_df <- map_dfr(colnames(coef_matrix), function(feature_name) {
  # Get the distribution of coefficients for the current feature
  feature_coeffs <- coef_matrix[, feature_name]

  # Calculate stability (frequency of non-zero coefficients)
  stability_freq <- sum(feature_coeffs != 0) / BOOTSTRAP_ITERATIONS

  # We are only interested in features with some minimum stability
  if (stability_freq >= MIN_FEATURE_STABILITY) {
    # Get the subset of non-zero coefficients to calculate CIs
    non_zero_coeffs <- feature_coeffs[feature_coeffs != 0]

    # Calculate the mean coefficient and the 95% CI from the bootstrap
    # distribution
    mean_coef <- mean(non_zero_coeffs)
    ci_coef <- quantile(
      x = non_zero_coeffs,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    # Exponentiate to calculate the odds ratio and its CI
    tibble(
      feature = feature_name,
      stability_freq = stability_freq,
      mean_coefficient = mean_coef,
      odds_ratio = exp(mean_coef),
      ci_low = exp(ci_coef[1]),
      ci_high = exp(ci_coef[2])
    )
  } else {
    NULL # Return nothing if the feature is not stable
  }
}) %>%
  arrange(desc(odds_ratio)) %>%
  # Is the odds ratio greater or less than one?
  mutate(
    direction = ifelse(
      odds_ratio > 1,
      yes = "greater",
      no = "lesser"
    )
  )

print(odds_ratio_df)
```

    ## # A tibble: 368 × 7
    ##    feature stability_freq mean_coefficient odds_ratio ci_low ci_high direction
    ##    <chr>            <dbl>            <dbl>      <dbl>  <dbl>   <dbl> <chr>    
    ##  1 vraS             0.998            1.54        4.65   2.17   11.8  greater  
    ##  2 liaF             1                1.44        4.22   2.12    9.77 greater  
    ##  3 rpsU             0.999            1.40        4.07   2.11    9.89 greater  
    ##  4 yycH             0.992            1.24        3.44   1.54    8.14 greater  
    ##  5 atl_1            0.965            1.04        2.84   1.25    7.21 greater  
    ##  6 pdhC_2           0.759            0.849       2.34   1.07    6.54 greater  
    ##  7 cdaA             0.588            0.845       2.33   1.08    8.35 greater  
    ##  8 lyrA             0.803            0.828       2.29   1.07    6.42 greater  
    ##  9 srrA             0.866            0.819       2.27   1.04    5.90 greater  
    ## 10 alaS             0.583            0.737       2.09   1.05    5.92 greater  
    ## # ℹ 358 more rows

``` r
# --- Process the odds ratio data for visualization and BLCA ---

# Exclude features in which the CI overlaps with 1 (i.e., not statistically
# significant), and exclude odds ratios with an absolute difference of more than
# 0.5
odds_ratio_processed <- odds_ratio_df %>%
  filter(!(ci_low <= 1 & ci_high >= 1)) %>%
  filter(abs(odds_ratio - 1) >= 0.5)

# Get features with odds ratio > 1 (associated with vancomycin)
vanc_features <- odds_ratio_processed %>%
  filter(direction == "greater") %>%
  pull(feature)

# Get features with odds ratio < 1 (associated with control)
control_features <- odds_ratio_processed %>%
  filter(direction == "lesser") %>%
  pull(feature)

# Processed features
processed_features <- odds_ratio_processed$feature

# Processed features with odds ratios greater than 1
processed_features_greater <- odds_ratio_processed %>%
  filter(odds_ratio > 1) %>%
  pull(feature)

# Processed features with odds ratios less than 1
processed_features_lesser <- odds_ratio_processed %>%
  filter(odds_ratio < 1) %>%
  pull(feature)

# Create a matrix of processed features for visualization and BLCA
processed_features_matrix <- mutation_matrix[, processed_features]
```

``` r
# --- Bar chart of mutation types (Figure 2A) ---

figure_2A <- variants %>%
  filter(is.na(notes) & !is.na(mutation_type)) %>%
  ggplot(
    aes(x = fct_infreq(mutation_type))
  ) +
  geom_bar() +
  xlab("Mutation type") +
  ylab("Count") +
  scale_x_discrete(limits = rev) +
  theme_cowplot() +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  coord_flip()


# --- Odds ratio plot (Figure 2B) ---

generate_OR_plot <- function(.x) {
  OR_filtered <- .x %>%
    filter(feature %in% processed_features)

  plot <- OR_filtered %>%
    ggplot(
      aes(
        x = fct_reorder(
          .f = feature,
          .x = odds_ratio,
          .desc = TRUE
        ),
        y = odds_ratio,
        color = direction
      )
    ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = ci_low,
        ymax = ci_high
      ),
      width = 0.2
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "#000000"
    ) +
    ylab("Average odds ratio") +
    xlab("Gene") +
    scale_color_manual(values = c("#847F9F", "#FAC9B8")) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(
        face = "italic",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "#EBEBEB"),
      legend.position = "none"
    )

  return(plot)
}

figure_2B <- generate_OR_plot(
  .x = odds_ratio_df
)


# --- Heatmap of mutations in selected features (Figure 2C) ---

generate_heatmap_plot <- function(.x) {
  # Prepare data for heatmap
  heatmap_df <- .x %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature") %>%
    pivot_longer(
      cols = -feature,
      names_to = "population",
      values_to = "mutation_presence"
    ) %>%
    mutate(
      direction = ifelse(
        test = feature %in% processed_features_greater,
        yes = "greater",
        no = "lesser"
      ),
      treatment_group = ifelse(
        test = grepl(
          pattern = "T",
          x = population
        ),
        yes = "vancomycin",
        no = "control"
      ),
      population = str_remove_all(
        string = population,
        pattern = "\\D"
      )
    )

  # Heatmap plot
  heatmap_plot <- heatmap_df %>%
    ggplot() +
    facet_grid(
      factor(
        x = treatment_group,
        levels = c("vancomycin", "control")
      ) ~ .,
      scales = "free_y",
      space = "free_y",
      labeller = as_labeller(
        c(
          "control" = "Control lines",
          "vancomycin" = "Vancomycin-adapted lines"
        )
      )
    ) +
    geom_tile(
      aes(
        y = fct_rev(
          factor(
            x = population,
            levels = 1:87
          )
        ),
        x = factor(
          x = feature,
          levels = colnames(processed_features_matrix)
        ),
        fill = mutation_presence
      ),
      data = filter(
        .data = heatmap_df,
        direction == "lesser"
      ),
      color = "#EBEBEB"
    ) +
    labs(
      y = "Population",
      x = "Gene",
    ) +
    scale_fill_gradient(low = "#FFFFFF", high = "#FAC9B8") +
    new_scale_fill() +
    geom_tile(
      aes(
        y = fct_rev(
          factor(
            x = population,
            levels = 1:87
          )
        ),
        x = factor(
          x = feature,
          levels = colnames(processed_features_matrix)
        ),
        fill = mutation_presence
      ),
      data = filter(
        .data = heatmap_df,
        direction == "greater"
      ),
      color = "#EBEBEB"
    ) +
    scale_x_discrete(limits = colnames(processed_features_matrix)) +
    scale_fill_gradient(low = "#FFFFFF", high = "#847F9F") +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_text(
        face = "italic",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )

  return(heatmap_plot)
}

figure_2C <- generate_heatmap_plot(.x = processed_features_matrix)

list(
  ... = figure_2A,
  ... = figure_2B,
  ... = figure_2C
)
```

    ## $...

![](S_aureus_evolution_files/figure-gfm/Plot%20of%20mutation%20heatmap-1.png)<!-- -->

    ## 
    ## $...

![](S_aureus_evolution_files/figure-gfm/Plot%20of%20mutation%20heatmap-2.png)<!-- -->

    ## 
    ## $...

![](S_aureus_evolution_files/figure-gfm/Plot%20of%20mutation%20heatmap-3.png)<!-- -->

We then perform a Bayesian latent class analysis (BLCA). This analysis
allows us to identify distinct mutation patterns associated with
vancomycin selection and control conditions. We perform the BLCA on the
processed mutation matrix, which includes only the features that were
selected by the weighted elastic-net logistic regression. The BLCA will
be run with a range of clusters to identify the optimal number based on
the AIC.

We use the results of the BLCA to construct a consensus matrix that
represents the frequency of co-occurrence of mutations across the
populations.

``` r
# --- Perform Bayesian Latent Class Analysis (BLCA) ---

perform_blca <- function(
    .x,
    group_label,
    iterations = BLCA_ITERATIONS,
    num_clusters = BLCA_CLUSTERS,
    num_cores = detectCores() - 1,
    output_dir_data = "data/cluster_results",
    output_dir_figs = "figures/cluster_plots") {


  # --- BLCA function ---

  blca_analysis <- function(data_matrix, group_label, num_clusters) {
    group_matrix <- data_matrix[
      grep(
        pattern = group_label,
        x = rownames(data_matrix)
      ),
    ]

    # We only include genes that have a mutation in at least 3 populations
    group_matrix <- group_matrix[, colSums(group_matrix) >= 3, drop = FALSE]

    # If the group matrix is too small to run BLCA, skip this iteration
    if (nrow(group_matrix) < 2 || ncol(group_matrix) < 1) {
      return(NULL)
    }

    # Create a local matrix by sampling 90% of the rows and columns
    local_matrix <- group_matrix[
      sample(
        x = seq_len(nrow(group_matrix)),
        size = ceiling(nrow(group_matrix) * 0.9),
        replace = FALSE
      ),
      sample(
        seq_len(ncol(group_matrix)),
        size = ceiling(ncol(group_matrix) * 0.9),
        replace = FALSE
      )
    ]

    tryCatch(
      { # Perform BLCA on the local matrix
        local_blca_res <- sapply(
          X = seq_len(num_clusters),
          simplify = FALSE,
          FUN = function(k) {
            blca.em(
              X = local_matrix,
              G = k,
              delta = 1.0,
              alpha = 0.5,
              beta = 0.5,
              restarts = 10,
              iter = 1000,
              start.vals = "across"
            )
          }
        )

        # Identify the best model based on AIC
        best_k <- which.max(
          sapply(
            X = local_blca_res,
            FUN = function(x) {
              x$AIC
            }
          )
        )

        best_res <- local_blca_res[[best_k]]

        # Cluster assignments
        if (best_k == 1) {
          cluster_assign <- setNames(
            rep(
              x = 1,
              times = nrow(local_matrix)
            ),
            rownames(local_matrix)
          )

          return(cluster_assign)

        } else if (!is.null(best_res$Z)) {
          clust_assign_uniq <- apply(
            X = best_res$Z,
            MARGIN = 1,
            FUN = which.max
          )

          # Need to match bitstring representations to observations
          temp <- local_matrix[, match(
            colnames(best_res$itemprob),
            table = colnames(local_matrix)
          )]

          # Match features just to be on the safe side
          bit_str <- apply(
            X = local_matrix,
            MARGIN = 1,
            FUN = function(s) {
              paste(s, collapse = "")
            }
          )

          cluster_assign <- setNames(
            clust_assign_uniq[match(bit_str, table = names(clust_assign_uniq))],
            rownames(local_matrix)
          )

          return(cluster_assign)

        } else {
          return(NULL)
        }
        return(cluster_assign)
      },
      warning = function(cond) {
        return("Warning")
      },
      error = function(cond) {
        return(cond)
      }
    )
  }

  # --- Parallelize the BLCA function ---

  parallelize_blca <- function(
    input = NULL,
    group_label = NULL,
    num_clusters = NULL
  ) {
    # Detect number of cores and register parallel backend

    cl <- makeCluster(num_cores)
    # Ensure required packages are available on workers
    clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(BayesLCA))
    })
    clusterExport(
      cl,
      varlist = c("blca_analysis", "input", "group_label", "num_clusters"),
      envir = environment()
    )

    blca_res_em <- parSapply(
      cl,
      1:iterations,
      simplify = FALSE,
      function(i) {
        temp <- blca_analysis(
          data_matrix = input,
          group_label = group_label,
          num_clusters = num_clusters
        )

        return(temp)
      }
    )

    stopCluster(cl)
    return(blca_res_em)
  }


  # --- Function to construct the consensus matrices ---

  construct_frequency_matrix <- function(
    data_matrix,
    blca_res_em,
    group_label
  ) {
    group_matrix <- data_matrix[
      grep(
        pattern = group_label, # Group label corresponds to control or treated
        x = rownames(data_matrix)
      ),
    ]

    sorted_ids <- sort(
      rownames(
        group_matrix
      )
    )

    pair_matrix <- Matrix(
      data = 0,
      ncol = length(sorted_ids),
      nrow = length(sorted_ids),
      dimnames = list(
        sorted_ids,
        sorted_ids
      )
    )

    count_matrix <- pair_matrix

    # Iterate over the LCA results to format the Consensus matrix.
    # The frequency is the number of co-occurrences divided by the number
    # of samples.

    for (i in seq_len(length(blca_res_em))) {
      local_res <- blca_res_em[[i]]

      # Skip entries that are NULL, errors (lists), warnings (characters),
      # or all NA
      if (is.null(local_res)) next
      if (!is.atomic(local_res)) next
      if (!is.numeric(local_res)) next
      if (length(local_res) == 0L || all(is.na(local_res))) next

      local_k <- suppressWarnings(max(local_res, na.rm = TRUE))

      if (!is.finite(local_k) || is.na(local_k)) {
        next
      }

      local_clust <- as.matrix(
        Matrix::Diagonal(n = local_k + 1)[local_res, ]
      )

      local_clust <- local_clust %*% t(local_clust)

      # Reorder
      idx <- match(
        x = sorted_ids,
        table = names(local_res)
      )

      local_clust <- local_clust[idx, idx]
      local_clust[is.na(local_clust)] <- 0
      pair_matrix <- pair_matrix + local_clust

      local_clust <- as.matrix(
        Matrix::Diagonal(n = 2)[rep(x = 1, length = length(local_res)), ]
      )

      local_clust <- local_clust %*% t(local_clust)

      idx <- match(
        sorted_ids,
        table = names(local_res)
      )

      local_clust <- local_clust[idx, idx]
      local_clust[is.na(local_clust)] <- 0
      count_matrix <- count_matrix + local_clust
    }

    freq_matrix <- pair_matrix # Initialize

    # Get indices where count_matrix is not zero
    valid_indices <- which(count_matrix > 0, arr.ind = TRUE)

    # Calculate frequency only for those indices
    freq_matrix[valid_indices] <- pair_matrix[valid_indices] / count_matrix[valid_indices]

    freq_matrix <- Matrix::forceSymmetric(freq_matrix)
    freq_matrix <- as.matrix(freq_matrix)

    return(freq_matrix)
  }


  # --- Hierarchical clustering and plotting function ---

  plot_clusters <- function(
    data_matrix,
    freq_matrix,
    group_label
  ) {
    group_matrix <- data_matrix[
      grep(
        pattern = group_label, # Group label corresponds to control or treated
        x = rownames(data_matrix)
      ),
    ]

    # Perform hierarchical clustering on the distance matrix
    dist_matrix <- as.dist(1.0 - freq_matrix)

    hc <- hclust(
      d = dist_matrix,
      method = "ward.D2"
    )

    output_tag <- ifelse(
      group_label == "^C",
      yes = "Control",
      no = "Vancomycin"
    )

    # Create directories if they don't exist
    dir.create(
      path = output_dir_data,
      showWarnings = FALSE,
      recursive = TRUE
    )

    dir.create(
      path = output_dir_figs,
      showWarnings = FALSE,
      recursive = TRUE
    )

    # Set colors for clusters
    cluster_colors <- setNames(
      RColorBrewer::brewer.pal(
        n = 7,
        name = "Paired"
      ),
      paste0("Cluster-", 1:7)
    )

    for (k in c(2:7)) {
      clust_assign <- cutree(
        tree = hc,
        k = k
      )

      clust_res <- data.frame(
        "ID" = names(clust_assign),
        "Cluster" = paste0("Cluster-", clust_assign)
      )

      write.csv(
        x = clust_res,
        file = file.path(
          output_dir_data,
          sprintf("ClusterResults_%s_K%s.csv", output_tag, k)
        ),
        quote = FALSE,
        row.names = FALSE
      )

      rowAnnot <- HeatmapAnnotation(
        Cluster = clust_res$Cluster[
          match(
            rownames(freq_matrix),
            table = clust_res$ID
          )
        ],
        col = list(Cluster = setNames(
          ggsci::pal_npg(palette = "nrc")(k),
          paste0("Cluster-", 1:k)
        )),
        which = "row",
        simple_anno_size = unit(3, "mm"),
        show_annotation_name = FALSE
      )

      col_fun_freq <- colorRamp2(
        breaks = seq(0, 1, 0.01),
        colors = colorRampPalette(
          colors = c("white", "#606B91")
        )(
          length(
            seq(0, 1, 0.01)
          )
        )
      )

      h_ae <- Heatmap(
        matrix = freq_matrix,
        col = col_fun_freq,
        cluster_rows = hc,
        cluster_columns = hc,
        show_column_names = TRUE,
        show_row_names = TRUE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        na_col = "white",
        heatmap_legend_param = list(
          title = "Consensus Frequency",
          title_position = "leftcenter-rot",
          legend_height = unit(4, "cm"),
          fontsize = 12
        ),
        show_heatmap_legend = TRUE,
        use_raster = FALSE,
        border = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(
          lwd = 0.5,
          col = "white"
        ),
        left_annotation = rowAnnot
      )

      col_fun <- colorRamp2(
        breaks = seq(
          from = 0,
          to = 1,
          by = 0.01
        ),
        colors = colorRampPalette(
          colors = c("white", "grey30")
        )(
          length(
            seq(
              from = 0,
              to = 1,
              by = 0.01
            )
          )
        )
      )

      h_ae_mut <- Heatmap(
        group_matrix[
          match(
            rownames(freq_matrix),
            table = rownames(group_matrix)
          ),
        ],
        col = col_fun,
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "manhattan",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        na_col = "white",
        show_heatmap_legend = FALSE,
        use_raster = FALSE,
        border = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(
          lwd = 0.5,
          col = "white"
        )
      )

      h_comb <- h_ae + h_ae_mut

      pdf(
        file = file.path(
          output_dir_figs,
          sprintf("CC_Genomic_%s_K%s.pdf", output_tag, k)
        ),
        width = 10,
        height = 6
      )

      draw(h_comb)
      dev.off()
    }
  }

  ## Function calls ##
  blca_res_em <- parallelize_blca(
    input = .x,
    group_label = group_label,
    num_clusters = num_clusters
  )

  freq_matrix <- construct_frequency_matrix(
    data_matrix = .x,
    blca_res_em = blca_res_em,
    group_label = group_label
  )

  plot_clusters(
    data_matrix = .x,
    freq_matrix = freq_matrix,
    group_label = group_label
  )

  return(NULL)
}


# Create group-specific matrices
vanc_blca_matrix <- mutation_matrix[, vanc_features]
control_blca_matrix <- mutation_matrix[, control_features]

# Run BLCA for the vancomycin group
perform_blca(
  .x = vanc_blca_matrix,
  iterations = BLCA_ITERATIONS,
  group_label = "^T" # Analyze only treated samples
)
```

    ## NULL

``` r
# Run BLCA for the control group
perform_blca(
  .x = control_blca_matrix,
  iterations = BLCA_ITERATIONS,
  group_label = "^C" # Analyze only control samples
)
```

    ## NULL

The BLCA showed that the 18 VISA populations clustered into two distinct
groups based on shared gene-level mutations. For instance, Cluster 1 was
uniquely characterized by mutations in *rpsU* while notably lacking
mutations in *yycH*. This mutual exclusivity was nonsignificant by
Fisher’s exact test (*P* = 0.0537). However, grouping of lines by the
functional level shows near complete mutual exclusivity of *rpsU* and
WalKR regulon mutations (Fisher’s exact test, *P* = 0.0049).

``` r
# --- Build contingency tables from processed_features_matrix
# (treated rows only) ---

# Helper: safe binary vector for a single gene feature
get_feature_vec <- function(mat, gene) {
  if (!gene %in% colnames(mat)) {
    warning(
      sprintf(
        "Feature '%s' not found in processed_features_matrix; treating as all zeros.", gene)
    )
    return(rep(0L, nrow(mat)))
  }
  as.integer(mat[, gene] > 0)
}

# Helper: binary vector if any gene in a set is present
get_any_vec <- function(mat, genes) {
  genes_present <- intersect(genes, colnames(mat))
  if (length(genes_present) == 0L) {
    warning("None of the specified genes were found; returning all zeros.")
    return(rep(0L, nrow(mat)))
  }
  as.integer(rowSums(mat[, genes_present, drop = FALSE] > 0) > 0)
}

# Helper: make a 2x2 table with consistent labeling
make_2x2 <- function(a_vec, b_vec, a_name, b_name) {
  stopifnot(length(a_vec) == length(b_vec))
  # Rows: a present/absent; Cols: b present/absent (to mirror prior labels)
  a_present_b_present <- sum(a_vec == 1 & b_vec == 1)
  a_present_b_absent  <- sum(a_vec == 1 & b_vec == 0)
  a_absent_b_present  <- sum(a_vec == 0 & b_vec == 1)
  a_absent_b_absent   <- sum(a_vec == 0 & b_vec == 0)
  out <- matrix(
    c(a_present_b_present, a_present_b_absent, a_absent_b_present, a_absent_b_absent),
    nrow = 2, byrow = TRUE
  )
  rownames(out) <- c(paste0(a_name, "_present"), paste0(a_name, "_absent"))
  colnames(out) <- c(paste0(b_name, "_present"), paste0(b_name, "_absent"))
  out
}

# Subset to vancomycin-treated populations (row names start with "T")
vanc_mat <- processed_features_matrix[grepl("^T", rownames(processed_features_matrix)), , drop = FALSE]

# 1) yycH vs rpsU
rpsU_vec <- get_feature_vec(vanc_mat, "rpsU")
yycH_vec <- get_feature_vec(vanc_mat, "yycH")
yycH_versus_rpsU <- make_2x2(a_vec = rpsU_vec, b_vec = yycH_vec, a_name = "rpsU", b_name = "yycH")
print(yycH_versus_rpsU)
```

    ##              yycH_present yycH_absent
    ## rpsU_present            0           6
    ## rpsU_absent             6           6

``` r
# 2) WalKR regulon (any of these genes) vs rpsU
walKR_regulon_genes <- c("walK", "atl_1", "yycH", "yycI")
regulon_vec <- get_any_vec(vanc_mat, walKR_regulon_genes)
regulon_versus_rpsU <- make_2x2(a_vec = rpsU_vec, b_vec = regulon_vec, a_name = "rpsU", b_name = "regulon")
print(regulon_versus_rpsU)
```

    ##              regulon_present regulon_absent
    ## rpsU_present               2              4
    ## rpsU_absent               12              0

``` r
list(
  yycH_vs_rpsU = fisher.test(yycH_versus_rpsU),
  regulon_vs_rpsU = fisher.test(regulon_versus_rpsU)
)
```

    ## $yycH_vs_rpsU
    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  yycH_versus_rpsU
    ## p-value = 0.05371
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.000000 1.377046
    ## sample estimates:
    ## odds ratio 
    ##          0 
    ## 
    ## 
    ## $regulon_vs_rpsU
    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  regulon_versus_rpsU
    ## p-value = 0.004902
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.000000 0.508626
    ## sample estimates:
    ## odds ratio 
    ##          0

------------------------------------------------------------------------

# Divergent genetic pathways under vancomycin exposure led to varied collateral responses

The following data report the MICs of 8 antibiotics – cefazolin,
clindamycin, daptomycin, gentamycin, meropenem, nafcillin,
trimethoprim-sulfamethoxazole, and vancomycin – on the MSSA ancestral
clone ATCC 29213 and evolved vancomycin-intermediate populations. We
quantified the collateral response of an evolved population as the
difference in its $\text{log}_2$-transformed MIC relative to the
ancestral clone. A population is collaterally resistant when its MIC is
higher than the ancestral MIC and collaterally sensitive when it is
lower.

We use these data to evaluate whether genomic differences that evolved
during vancomycin adaptation shape drug tradeoffs.

``` r
collateral_MICs <- read_csv(
  file = input_collateral_MICs,
  show_col_types = FALSE
)
```

``` r
# --- Function that calculates collateral response values ---

calculate_collateral_responses <- function(.x) {
  # This function takes a data frame of MICs and calculates drug response
  # values for each evolved population relative to the ancestral population.
  ancestor_frame <- .x %>%
    filter(population == "Ancestor")

  evolved_frame <- .x %>%
    filter(population != "Ancestor")

  CR_values <- c()

  # This code chunk iterates through the unique paired_IDs in the
  # evolved_frame and filters the ancestor_frame and evolved_frame for the
  # paired_ID. The MIC values for the ancestor and evolved populations are
  # then used to calculate the collateral response value for each evolved
  # population. The collateral response values are then appended to the
  # CR_values vector.
  for (i in unique(evolved_frame$paired_ID)) {
    ancestor_value <- ancestor_frame %>%
      filter(paired_ID == i) %>%
      pull(MIC)

    evolved_subframe <- evolved_frame %>%
      filter(paired_ID == i)

    for (j in 1:nrow(evolved_subframe)) {
      evolved_value <- evolved_subframe[j, ] %>%
        pull(MIC)

      CR_values <- c(CR_values, log2(evolved_value / ancestor_value))
    }
  }

  # Creates a tibble of collateral response values
  CR_col <- tibble(CR = CR_values)
  CR_df <- bind_cols(
    ... = evolved_frame,
    ... = CR_col
  )

  # Remove paired_ID column
  CR_df <- CR_df %>%
    dplyr::select(-paired_ID)

  return(CR_df)
}


# --- Function that performs summary statistics on collateral responses ---

perform_summary_stats <- function(.x) {
  summary_stats <- .x %>%
    group_by(antibiotic, population) %>%
    summarize(
      mean_MIC = mean(MIC),
      mean_CR = mean(CR)
    )

  # Create a new column called "state" in the summary_stats data frame.
  # Populate this column with "-1", "0", or " 1" if the mean collateral
  # response value is less than, equal to, or greater than 0.
  summary_stats$state <- ifelse(
    test = summary_stats$mean_CR < 0,
    yes = "-1",
    no = ifelse(
      test = summary_stats$mean_CR == 0,
      yes = "0",
      no = "1"
    )
  )

  return(summary_stats)
}


# Calculates collateral responses
collateral_responses <- calculate_collateral_responses(.x = collateral_MICs)

# Calculates summary statistics
CR_summary_stats <- perform_summary_stats(.x = collateral_responses)

list(
  collateral_responses = collateral_responses,
  collateral_responses_summary = CR_summary_stats
)
```

    ## $collateral_responses
    ## # A tibble: 432 × 5
    ##    population replicate antibiotic   MIC     CR
    ##    <chr>          <dbl> <chr>      <dbl>  <dbl>
    ##  1 1                  1 CFZ         0.25 -0.585
    ##  2 1                  2 CFZ         0.25 -0.585
    ##  3 1                  3 CFZ         0.25 -0.585
    ##  4 2                  1 CFZ         0.5   0.415
    ##  5 2                  2 CFZ         0.5   0.415
    ##  6 2                  3 CFZ         0.5   0.415
    ##  7 3                  1 CFZ         0.75  1    
    ##  8 3                  2 CFZ         0.75  1    
    ##  9 3                  3 CFZ         0.75  1    
    ## 10 4                  1 CFZ         0.25 -0.585
    ## # ℹ 422 more rows
    ## 
    ## $collateral_responses_summary
    ## # A tibble: 144 × 5
    ## # Groups:   antibiotic [8]
    ##    antibiotic population mean_MIC mean_CR state
    ##    <chr>      <chr>         <dbl>   <dbl> <chr>
    ##  1 CFZ        1             0.25   -0.585 -1   
    ##  2 CFZ        10            0.5     0.138 1    
    ##  3 CFZ        11            0.5     0     0    
    ##  4 CFZ        12            0.583   0.333 1    
    ##  5 CFZ        13            0.75    1     1    
    ##  6 CFZ        14            0.5     0.415 1    
    ##  7 CFZ        15            0.458   0.277 1    
    ##  8 CFZ        16            0.25   -0.585 -1   
    ##  9 CFZ        17            0.375  -0.277 -1   
    ## 10 CFZ        18            0.5     0     0    
    ## # ℹ 134 more rows

``` r
figure_3B <- CR_summary_stats %>%
  ggplot(
    aes(
      y = population,
      x = antibiotic,
      fill = mean_CR
    )
  ) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "#1a4891",
    mid = "white",
    high = "#a01415",
    midpoint = 0
  ) +
  scale_y_discrete(
    limits = rev(
      factor(
        seq(1, 18, 1)
      )
    )
  ) +
  scale_x_discrete(
    limits =
      factor(
        c("DAP", "VAN", "CFZ", "SXT", "CLI", "MEM", "NAF", "GEN")
      ),
    position = "top"
  ) +
  scale_y_discrete(
    limits =
      factor(
        rev(
          c(4, 5, 6, 7, 8, 11, 14, 15, 16, 17, 18, 1, 2, 3, 10, 12, 13, 9)
        )
      )
  ) +
  ylab("Population") +
  xlab("Antibiotic") +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.ticks = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
# --- Function to plot collateral responses (Figure 3C) ---

generate_CR_plot <- function(.x) {
  # Factor antibiotics by order of decreasing median collateral response
  antibiotic_factor <- .x %>%
    group_by(antibiotic) %>%
    summarize(med = median(mean_CR)) %>%
    arrange(
      desc(med)
    ) %>%
    pull(antibiotic)

  .x$antibiotic <- factor(
    x = .x$antibiotic,
    levels = as.factor(antibiotic_factor)
  )

  CR_boxplot <- .x %>%
    ggplot(
      aes(
        x = antibiotic,
        y = mean_CR,
        fill = antibiotic
      )
    ) +
    geom_boxplot(outlier.shape = 1) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "#000000"
    ) +
    xlab("Antibiotic") +
    ylab(
      expression(
        "Log"[2] ~ "-" ~ "relative MIC"
      )
    ) +
    scale_fill_manual(
      values = c(
        "#db8f81", "#df998b", "#fae7e3", "#FFFFFF",
        "#f2f3f8", "#ececf5", "#dbddec", "#d1d4e7"
      )
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white")
    )

  return(CR_boxplot)
}


figure_3C <- generate_CR_plot(
  .x = CR_summary_stats
)

list(
  figure_3B = figure_3B,
  figure_3C = figure_3C
)
```

    ## $figure_3B

![](S_aureus_evolution_files/figure-gfm/Plot%20collateral%20response%20values-1.png)<!-- -->

    ## 
    ## $figure_3C

![](S_aureus_evolution_files/figure-gfm/Plot%20collateral%20response%20values-2.png)<!-- -->

We used Mann-Whitney U tests to compare the collateral response values
of the evolved populations to the ancestral population for each
antibiotic.

``` r
# --- Function that performs Mann-Whitney U tests on collateral response
# values ---

perform_MW_test <- function(.data, ab) {
  ancestor_df <- .data %>%
    filter(population == "Ancestor" & antibiotic == ab) %>%
    dplyr::select(
      ... = -replicate,
      ... = -paired_ID
    )

  evolved_df <- .data %>%
    filter(population != "Ancestor" & antibiotic == ab) %>%
    group_by(population) %>%
    # Calculate the mean MIC for each evolved population
    summarize(MIC = mean(MIC)) %>%
    # Relabel all entries in the "population" column to "Evolved"
    mutate(
      population = ifelse(
        population == "Ancestor",
        yes = "Ancestor",
        no = "Evolved"
      )
    )

  # Bind the antibiotic_df and ancestor_df objects
  bound_df <- bind_rows(
    ... = evolved_df,
    ... = ancestor_df
  )

  # Perform Mann Whitney U test comparing the MICs of the ancestor and
  # evolved populations
  MW_test <- wilcox.test(
    formula = MIC ~ population,
    data = bound_df,
    exact = FALSE,
    correct = FALSE,
    alternative = "two.sided"
  )

  return(MW_test)
}

mw_stats <- map(
  unique(collateral_MICs$antibiotic),
  ~ perform_MW_test(
    .data = collateral_MICs,
    ab = .x
  )
)

names(mw_stats) <- unique(collateral_MICs$antibiotic)

list(mw_stats)
```

    ## [[1]]
    ## [[1]]$CFZ
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 55, p-value = 0.3321
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$CLI
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 108, p-value = 0.01754
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$DAP
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 4, p-value = 0.0001314
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$GEN
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 97.5, p-value = 0.1506
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$MEM
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 84, p-value = 0.4946
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$NAF
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 96.5, p-value = 0.156
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$SXT
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 44, p-value = 0.04432
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## [[1]]$VAN
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  MIC by population
    ## W = 0, p-value = 3.849e-05
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# Extract p-values from the mw_stats object and put into a dataframe
mw_pvals <- c()
for (i in 1:length(mw_stats)) {
  mw_pvals <- rbind(
    mw_pvals,
    data.frame(
      antibiotic = names(mw_stats)[i],
      p_value = mw_stats[[i]]$p.value
    )
  )
}

mw_p_adjusted <- p.adjust(
  p = mw_pvals$p_value,
  method = "BH"
)

names(mw_p_adjusted) <- mw_pvals$antibiotic

mw_p_adjusted
```

    ##          CFZ          CLI          DAP          GEN          MEM          NAF 
    ## 0.3794910158 0.0467799105 0.0005256211 0.2080483531 0.4946242449 0.2080483531 
    ##          SXT          VAN 
    ## 0.0886339265 0.0003079510

We can combine the collateral response and genomic data obtained from
each population to evaluate whether genomic differences that evolved
during vancomycin adaptation shape drug tradeoffs.

``` r
# --- Function to perform biserial correlation analysis ---

biserial_correlation <- function(.data, lines, ab) {
  if (length(lines) == 12) {
    .data <- .data %>%
      filter(!population %in% c(1, 2))
  }

  .data$mutated <- ifelse(
    .data$population %in% lines,
    "1",
    "0"
  )

  ab_frame <- .data %>%
    filter(antibiotic == ab)

  mutated <- ab_frame %>%
    pull(mutated)

  mean_CR <- ab_frame %>%
    pull(mean_CR)

  results <- cor.test(
    x = mean_CR,
    y = as.numeric(mutated)
  )

  return(results)
}

# Helper function
perform_biserial_correlation <- function(.data) {
  antibiotic_vector <- c("CFZ", "MEM", "NAF", "GEN")

  regulon_lines <- c(4, 5, 6, 7, 8, 9, 11, 14, 15, 16, 17, 18)
  yycHI_lines <- c(4, 5, 6, 7, 9, 16, 18)
  rpsU_lines <- c(3, 10, 12, 13)

  regulon_results <- map(
    .x = antibiotic_vector,
    ~ biserial_correlation(
      .data = .data,
      lines = regulon_lines,
      ab = .x
    )
  )

  yycHI_results <- map(
    .x = antibiotic_vector,
    ~ biserial_correlation(
      .data = .data,
      lines = yycHI_lines,
      ab = .x
    )
  )

  rpsU_results <- map(
    .x = antibiotic_vector,
    ~ biserial_correlation(
      .data = .data,
      lines = rpsU_lines,
      ab = .x
    )
  )

  names(regulon_results) <- antibiotic_vector
  names(yycHI_results) <- antibiotic_vector
  names(rpsU_results) <- antibiotic_vector

  return(
    list(
      regulon_results = regulon_results,
      yycHI_results = yycHI_results,
      rpsU_results = rpsU_results
    )
  )
}

biserial_results <- perform_biserial_correlation(.data = CR_summary_stats)


regulon_p_adjusted <- p.adjust(
  p = map_dbl(biserial_results$regulon_results, ~ .x$p.value),
  method = "BH"
)

rpsU_p_adjusted <- p.adjust(
  p = map_dbl(biserial_results$rpsU_results, ~ .x$p.value),
  method = "BH"
)

list(
  biserial_results = biserial_results,
  regulon_p_adjusted = regulon_p_adjusted,
  yycH_p = biserial_results$yycH_results$NAF$p.value,
  rpsU_p_adjusted = rpsU_p_adjusted
)
```

    ## $biserial_results
    ## $biserial_results$regulon_results
    ## $biserial_results$regulon_results$CFZ
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -2.1789, df = 14, p-value = 0.04692
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.79949543 -0.01001952
    ## sample estimates:
    ##        cor 
    ## -0.5032255 
    ## 
    ## 
    ## $biserial_results$regulon_results$MEM
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -1.9766, df = 14, p-value = 0.06812
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.78178927  0.03721915
    ## sample estimates:
    ##        cor 
    ## -0.4671041 
    ## 
    ## 
    ## $biserial_results$regulon_results$NAF
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -3.3833, df = 14, p-value = 0.004458
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8753694 -0.2621456
    ## sample estimates:
    ##       cor 
    ## -0.670696 
    ## 
    ## 
    ## $biserial_results$regulon_results$GEN
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = 1.1847, df = 14, p-value = 0.2559
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2279566  0.6937553
    ## sample estimates:
    ##       cor 
    ## 0.3018586 
    ## 
    ## 
    ## 
    ## $biserial_results$yycHI_results
    ## $biserial_results$yycHI_results$CFZ
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -2.4982, df = 16, p-value = 0.02376
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.79899078 -0.08349926
    ## sample estimates:
    ##        cor 
    ## -0.5297191 
    ## 
    ## 
    ## $biserial_results$yycHI_results$MEM
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -0.97185, df = 16, p-value = 0.3456
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.6331719  0.2593655
    ## sample estimates:
    ##       cor 
    ## -0.236093 
    ## 
    ## 
    ## $biserial_results$yycHI_results$NAF
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -2.2719, df = 16, p-value = 0.03724
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.78073109 -0.03510472
    ## sample estimates:
    ##        cor 
    ## -0.4938804 
    ## 
    ## 
    ## $biserial_results$yycHI_results$GEN
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = 1.6175, df = 16, p-value = 0.1253
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1115013  0.7163729
    ## sample estimates:
    ##      cor 
    ## 0.374884 
    ## 
    ## 
    ## 
    ## $biserial_results$rpsU_results
    ## $biserial_results$rpsU_results$CFZ
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = 2.2404, df = 16, p-value = 0.03961
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.02823668 0.77803225
    ## sample estimates:
    ##       cor 
    ## 0.4886648 
    ## 
    ## 
    ## $biserial_results$rpsU_results$MEM
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = 2.1195, df = 16, p-value = 0.05004
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.00170504 0.76733989
    ## sample estimates:
    ##       cor 
    ## 0.4682024 
    ## 
    ## 
    ## $biserial_results$rpsU_results$NAF
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = 2.9472, df = 16, p-value = 0.009465
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.1746891 0.8301505
    ## sample estimates:
    ##      cor 
    ## 0.593181 
    ## 
    ## 
    ## $biserial_results$rpsU_results$GEN
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  mean_CR and as.numeric(mutated)
    ## t = -0.86785, df = 16, p-value = 0.3983
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.6177480  0.2828396
    ## sample estimates:
    ##        cor 
    ## -0.2120288 
    ## 
    ## 
    ## 
    ## 
    ## $regulon_p_adjusted
    ##        CFZ        MEM        NAF        GEN 
    ## 0.09083165 0.09083165 0.01783375 0.25585266 
    ## 
    ## $yycH_p
    ## NULL
    ## 
    ## $rpsU_p_adjusted
    ##        CFZ        MEM        NAF        GEN 
    ## 0.06672216 0.06672216 0.03785965 0.39830806

``` r
figure_3D <- CR_summary_stats %>%
  filter(antibiotic %in% c("CFZ", "MEM", "NAF", "GEN")) %>%
  # Remove populations that have both regulon and rpsU mutations
  filter(!population %in% c(1, 2)) %>%
  mutate(
    regulon = ifelse(
      population %in% c(4, 5, 6, 7, 8, 9, 11, 14, 15, 16, 17, 18),
      "1",
      "0"
    )
  ) %>%
  ggplot(
    aes(
      x = antibiotic,
      y = mean_CR,
      fill = regulon
    )
  ) +
  geom_boxplot(outlier.shape = 1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#000000"
  ) +
  xlab("Antibiotic") +
  ylab(
    expression(
      "Log"[2] ~ "-" ~ "relative MIC"
    )
  ) +
  scale_fill_manual(
    labels = c("rpsU mutations", "Regulon mutations"),
    values = c("#FFFFFF", "#C2C2C2")
  ) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "white"),
    legend.title = element_blank()
  )

figure_3D
```

![](S_aureus_evolution_files/figure-gfm/Regulon%20plot-1.png)<!-- -->

``` r
figure_3E <- CR_summary_stats %>%
  filter(antibiotic == "NAF") %>%
  mutate(
    mutation = ifelse(
      population %in% c(4, 5, 6, 7, 9, 16, 18),
      "1",
      "0"
    )
  ) %>%
  ggplot(
    aes(
      x = antibiotic,
      y = mean_CR,
      fill = mutation
    )
  ) +
  geom_boxplot(outlier.shape = 1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#000000"
  ) +
  xlab("Antibiotic") +
  ylab(
    expression(
      "Log"[2] ~ "-" ~ "relative MIC"
    )
  ) +
  scale_fill_manual(
    labels = c("No yycH mutations", "yycH mutations"),
    values = c("#FFFFFF", "#C2C2C2")
  ) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    axis.title.x = element_blank()
  )

figure_3E
```

![](S_aureus_evolution_files/figure-gfm/yycH%20plot-1.png)<!-- -->

------------------------------------------------------------------------

# Using Collateral Response Scores to inform antibiotic treatment

Now, we first estimate the bootstrap distribution of collateral response
scores (CRS) for each antibiotic and mutational pathway. The CSS metric
evaluates the net collateral effect of antibiotic exposure by
integrating both the magnitude and direction of changes in MIC,
providing a robust metric for evaluating the likelihood of collateral
sensitivity or resistance. This score is amenable to statistical
investigation and provides a standardized way of reporting this
information for clinical decision-making. The CSS is defined as,

$$\text{CSS} = \frac{\frac{1}{N} \sum_{i=1}^{N} \log_2\left
(\frac{\mathrm{MIC}_{i,\text{evolved}}}{\mathrm{MIC}_
{i,\text{ancestor}}}\right)}{\max_{1 \leq i \leq N} \left|\log_2\left
(\frac{\mathrm{MIC}_{i,\text{evolved}}}{\mathrm{MIC}_
{i,\text{ancestor}}}\right)\right|}$$

where *N* is the number of replicate populations, *i*.

``` r
# Function that computes bootstrap collateral response scores (CRS) and
# percentile-based confidence intervals for each drug and mutational pathway
# and plots these distributions.

perform_bootstrap_and_plot_css <- function(
  .x,
  groups,
  iterations = CRS_ITERATIONS
) {
  # Function to perform bootstrap analysis
  perform_bootstrap <- function(data, groups, iterations) {
    # CSS = (mean collateral response / max |collateral response|)
    css_stat <- function(data, indices) {
      d <- data[indices]
      mean_val <- mean(d)
      max_val <- max(abs(d))
      css <- mean_val / max_val
      return(css)
    }

    resampling <- function(df, R) {
      boot_obj <- boot(
        data = df$mean_CR,
        statistic = css_stat,
        R = R
      )

      # Extract raw bootstrap replicates
      raw_boot <- data.frame(CSS = boot_obj$t[, 1])
      raw_boot$replicate <- seq_len(nrow(raw_boot))

      for (g in groups) {
        raw_boot[[g]] <- unique(df[[g]])
      }

      # Compute percentile-based confidence intervals
      ci_obj <- boot.ci(
        boot.out = boot_obj,
        type = "perc"
      )

      lower_perc <- ci_obj$percent[4]
      upper_perc <- ci_obj$percent[5]

      # Compute point estimate
      point_est <- css_stat(
        data = df$mean_CR,
        indices = 1:length(df$mean_CR)
      )

      raw_boot$lower_perc <- lower_perc
      raw_boot$upper_perc <- upper_perc
      raw_boot$point_est <- point_est

      return(raw_boot)
    }

    results <- data %>%
      group_by(across(all_of(groups))) %>%
      do(
        resampling(
          df = .,
          R = iterations
        )
      ) %>%
      ungroup()

    return(results)
  }

  # Function to create three ridgeline plots given the original data.
  # It computes bootstrap results for:
  # 1. Grouping by antibiotic only,
  # 2. Grouping by cluster only,
  # 3. Grouping by both antibiotic and cluster.

  generate_plots <- function(data, iterations) {
    # Bootstrap results for different grouping schemes.
    boot_res_ab <- perform_bootstrap(
      data = data,
      groups = c("antibiotic"),
      iterations = iterations
    )

    boot_res_cluster <- perform_bootstrap(
      data = data,
      groups = c("cluster"),
      iterations = iterations
    )

    boot_res_both <- perform_bootstrap(
      data = data,
      groups = c("antibiotic", "cluster"),
      iterations = iterations
    )

    # Create a numeric representation and offset only for segments/points.
    boot_res_both <- boot_res_both %>%
      mutate(
        antibiotic_numeric = as.numeric(
          factor(
            x = antibiotic,
            levels = unique(antibiotic)
          )
        ),
        y_offset = antibiotic_numeric + if_else(as.factor(cluster) == "0", -0.2, -0.1)
      )

    # Summarize data for segments and points.
    summary_data <- boot_res_both %>%
      group_by(antibiotic, cluster) %>%
      summarize(
        lower = unique(lower_perc),
        upper = unique(upper_perc),
        point_est = unique(point_est),
        y_offset = unique(y_offset),
        .groups = "drop"
      )

    ## Plot 1: By antibiotic only.
    plot_ab <- ggplot(
      data = boot_res_ab,
      aes(
        x = CSS,
        y = fct_reorder(antibiotic, point_est)
      )
    ) +
      geom_density_ridges(
        scale = 1.2,
        fill = "#c2c2c2",
        color = "black"
      ) +
      geom_segment(
        data = boot_res_ab %>%
          group_by(antibiotic) %>%
          summarize(
            lower = unique(lower_perc),
            upper = unique(upper_perc),
            .groups = "drop"
          ),
        aes(
          x = lower,
          xend = upper,
          y = antibiotic,
          yend = antibiotic
        ),
        color = "#7A7A7A",
        linewidth = 1
      ) +
      geom_point(
        data = boot_res_ab %>%
          group_by(antibiotic) %>%
          summarize(
            mean_css = unique(point_est),
            .groups = "drop"
          ),
        aes(
          x = mean_css,
          y = antibiotic
        ),
        color = "#7A7A7A",
        size = 3
      ) +
      scale_x_continuous(
        limits = c(-1, 1),
        breaks = seq(-1, 1, 0.25)
      ) +
      labs(
        x = "Collateral Response Score (CSS)",
        y = "Antibiotic"
      ) +
      theme_cowplot() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 12)
      )

    ## Plot 2: By cluster only.
    plot_cluster <- ggplot(
      data = boot_res_cluster,
      aes(
        x = CSS,
        y = fct_reorder(
          as.factor(cluster),
          point_est
        ),
        fill = as.factor(cluster)
      )
    ) +
      geom_density_ridges(
        alpha = 0.5,
        scale = 1.2,
        color = "black"
      ) +
      geom_segment(
        data = boot_res_cluster %>%
          group_by(cluster) %>%
          summarize(
            lower = unique(lower_perc),
            upper = unique(upper_perc),
            .groups = "drop"
          ),
        aes(
          x = lower,
          xend = upper,
          y = cluster,
          yend = cluster,
          color = as.factor(cluster)
        ),
        linewidth = 1
      ) +
      geom_point(
        data = boot_res_cluster %>%
          group_by(cluster) %>%
          summarize(
            mean_css = unique(point_est),
            .groups = "drop"
          ),
        aes(
          x = mean_css,
          y = cluster,
          color = as.factor(cluster)
        ),
        size = 3
      ) +
      labs(
        x = "Collateral Response Score (CRS)",
        y = "Cluster"
      ) +
      scale_x_continuous(
        limits = c(-1, 1),
        breaks = seq(-1, 1, 0.25)
      ) +
      scale_y_discrete(labels = rev(c("rpsU mutations", "Regulon mutations"))) +
      scale_color_manual(
        values = c("#83B9D8", "#A2B9B5")
      ) +
      scale_fill_manual(
        values = c("#83B9D8", "#A2B9B5")
      ) +
      theme_cowplot() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 12)
      )

    ## Plot 3: By antibiotic and cluster.
    plot_both <- ggplot() +
      # Density ridges centered on the antibiotic factor.
      geom_density_ridges(
        data = boot_res_both,
        aes(
          x = CSS,
          y = antibiotic, # use factor for centering
          fill = as.factor(cluster)
        ),
        alpha = 0.5,
        scale = 1.2,
        color = "black"
      ) +
      # Segments offset using y_offset.
      geom_segment(
        data = summary_data,
        aes(
          x = lower,
          xend = upper,
          y = y_offset,
          yend = y_offset,
          color = as.factor(cluster)
        ),
        linewidth = 1
      ) +
      # Points offset using y_offset.
      geom_point(
        data = summary_data,
        aes(
          x = point_est,
          y = y_offset,
          color = as.factor(cluster)
        ),
        size = 3
      ) +
      scale_fill_manual(
        values = c("#83B9D8", "#A2B9B5"),
        labels = c("rpsU mutations", "Regulon mutations")
      ) +
      scale_color_manual(
        values = c("#83B9D8", "#A2B9B5")
      ) +
      scale_x_continuous(
        limits = c(-1, 1),
        breaks = seq(-1, 1, 0.25)
      ) +
      labs(
        x = "Collateral Response Score (CSS)",
        y = "Antibiotic",
        fill = "Cluster"
      ) +
      theme_cowplot() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12)
      )

    return(
      results_list <- list(
        antibiotic_results = boot_res_ab,
        cluster_results = boot_res_cluster,
        both_results = boot_res_both,
        antibiotic_plot = plot_ab,
        cluster_plot = plot_cluster,
        both_plot = plot_both
      )
    )
  }

  # WalKR regulon
  cluster <- c(4, 5, 6, 7, 8, 9, 11, 14, 15, 16, 17, 18)

  .x <- .x %>%
    filter(!population %in% c(1, 2)) # Remove populations 1 and 2

  .x$cluster <- ifelse(
    .x$population %in% cluster,
    "1", # WalKR regulon
    "0" # rpsU mutations
  )

  plots <- generate_plots(
    data = .x,
    iterations = iterations
  )

  return(plots)
}

css_bootstrap_results <- perform_bootstrap_and_plot_css(.x = CR_summary_stats)

css_bootstrap_results
```

    ## $antibiotic_results
    ## # A tibble: 800 × 6
    ##        CSS replicate antibiotic lower_perc upper_perc point_est
    ##      <dbl>     <int> <chr>           <dbl>      <dbl>     <dbl>
    ##  1  0.206          1 CFZ            -0.172      0.387    0.0357
    ##  2  0.178          2 CFZ            -0.172      0.387    0.0357
    ##  3  0.0511         3 CFZ            -0.172      0.387    0.0357
    ##  4  0.129          4 CFZ            -0.172      0.387    0.0357
    ##  5 -0.134          5 CFZ            -0.172      0.387    0.0357
    ##  6  0.299          6 CFZ            -0.172      0.387    0.0357
    ##  7 -0.0280         7 CFZ            -0.172      0.387    0.0357
    ##  8  0.0389         8 CFZ            -0.172      0.387    0.0357
    ##  9  0.0806         9 CFZ            -0.172      0.387    0.0357
    ## 10  0.0764        10 CFZ            -0.172      0.387    0.0357
    ## # ℹ 790 more rows
    ## 
    ## $cluster_results
    ## # A tibble: 200 × 6
    ##      CSS replicate cluster lower_perc upper_perc point_est
    ##    <dbl>     <int> <chr>        <dbl>      <dbl>     <dbl>
    ##  1 0.198         1 0           0.0405      0.333     0.183
    ##  2 0.248         2 0           0.0405      0.333     0.183
    ##  3 0.194         3 0           0.0405      0.333     0.183
    ##  4 0.230         4 0           0.0405      0.333     0.183
    ##  5 0.303         5 0           0.0405      0.333     0.183
    ##  6 0.142         6 0           0.0405      0.333     0.183
    ##  7 0.209         7 0           0.0405      0.333     0.183
    ##  8 0.250         8 0           0.0405      0.333     0.183
    ##  9 0.314         9 0           0.0405      0.333     0.183
    ## 10 0.191        10 0           0.0405      0.333     0.183
    ## # ℹ 190 more rows
    ## 
    ## $both_results
    ## # A tibble: 1,600 × 9
    ##      CSS replicate antibiotic cluster lower_perc upper_perc point_est
    ##    <dbl>     <int> <chr>      <chr>        <dbl>      <dbl>     <dbl>
    ##  1 0.618         1 CFZ        0            0.403          1     0.618
    ##  2 0.618         2 CFZ        0            0.403          1     0.618
    ##  3 0.569         3 CFZ        0            0.403          1     0.618
    ##  4 1             4 CFZ        0            0.403          1     0.618
    ##  5 0.618         5 CFZ        0            0.403          1     0.618
    ##  6 0.785         6 CFZ        0            0.403          1     0.618
    ##  7 0.403         7 CFZ        0            0.403          1     0.618
    ##  8 0.569         8 CFZ        0            0.403          1     0.618
    ##  9 0.785         9 CFZ        0            0.403          1     0.618
    ## 10 0.667        10 CFZ        0            0.403          1     0.618
    ## # ℹ 1,590 more rows
    ## # ℹ 2 more variables: antibiotic_numeric <dbl>, y_offset <dbl>
    ## 
    ## $antibiotic_plot

    ## Picking joint bandwidth of 0.0354

![](S_aureus_evolution_files/figure-gfm/Compute%20CSS%20distributions-1.png)<!-- -->

    ## 
    ## $cluster_plot

    ## Picking joint bandwidth of 0.0179

![](S_aureus_evolution_files/figure-gfm/Compute%20CSS%20distributions-2.png)<!-- -->

    ## 
    ## $both_plot

    ## Picking joint bandwidth of 0.05

    ## Warning: Removed 45 rows containing non-finite outside the scale range
    ## (`stat_density_ridges()`).

![](S_aureus_evolution_files/figure-gfm/Compute%20CSS%20distributions-3.png)<!-- -->
