### Functions associated with the Card and Crozier et al. (2025) R Notebook
### Authors: Kyle Card and Arda Durmaz


## -----------------------------------------------------------------------------
## Function that plots vancomycin MIC data with the experimental schema
## (Figure 1)
## -----------------------------------------------------------------------------

generate_figure_1 <- function(.x) {
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

  # Combine the experimental schematic illustration and vancomycin plot
  # Imports PDF experimental schema image into R and adds it above the
  # vancomycin_MICs_plot
  schema <- image_read_pdf(
    "Figures/experimental_schema.pdf",
    density = 300
  )

  schema_grob <- rasterGrob(
    as.raster(schema),
    interpolate = TRUE
  )

  combined_plot <- plot_grid(
    schema_grob,
    vanc_plot,
    ncol = 1,
    rel_heights = c(0.3, 0.7),
    labels = c("A", "B"),
    label_size = 16
  )

  combined_plot <- combined_plot +
    theme(plot.background = element_rect(fill = "#FFFFFF"))

  return(combined_plot)
}


## -----------------------------------------------------------------------------
## Function that generates a binary matrix of qualifying gene-level mutations.
## -----------------------------------------------------------------------------

generate_binary_matrix <- function(.x) {
  # Creates a vector of the unique sample designations in the
  # qualifying_mutations object
  populations <- unique(.x$population)

  # Sorts the genes with qualifying_mutations by number of hits and
  # creates a vector
  genes <- .x %>%
    group_by(gene) %>%
    summarize(hits = n()) %>%
    arrange(
      desc(hits)
    ) %>%
    pull(gene)

  # Creates a binary matrix with the number of rows equal to the number of
  # samples, and the number of columns equal to the number of genes in the
  # mutations object and populates the matrix with 0s
  binary_mat <- matrix(
    data = 0,
    nrow = length(populations),
    ncol = length(genes)
  )

  rownames(binary_mat) <- populations
  colnames(binary_mat) <- genes

  # Use the qualifying_mutations object to populate the binary_mat object.
  # This code chunk iterates through the qualifying mutations data frame
  # and populates the binary_mat object with a 1 if a mutation is present in a
  # given gene / lineage combination
  for (i in seq_len(nrow(.x))) {
    binary_mat[.x$population[i], .x$gene[i]] <- 1
  }

  return(binary_mat)
}


## -----------------------------------------------------------------------------
## Function that performs similarity analysis and a permutation test on the
## qualifying gene-level mutations in the vancomycin-treated and control lines.
## -----------------------------------------------------------------------------

perform_genomic_analysis <- function(.x, num_perm = 10000) {
  calculate_similarity <- function(.x, permutation = FALSE) {
    # Sample the lines without replacement (shuffles population labels)
    if (permutation == TRUE) {
      shuffled_row_names <- sample(row.names(.x))
      row.names(.x) <- shuffled_row_names
    }

    # Compute Dice's Similarity Coefficient for all possible replicate pairs.
    # Populates a matrix with these values
    pwise_simil_matrix <- as.matrix(
      simil(
        x = .x,
        method = "Dice",
        by_rows = TRUE,
        upper = TRUE,
        diag = TRUE
      )
    )

    # All values above (and including) the matrix diagonal are converted to
    # NA to ease downstream wrangling
    diag(pwise_simil_matrix) <- NA
    pwise_simil_matrix[upper.tri(pwise_simil_matrix)] <- NA

    # Converted matrix to a data frame for analysis. NA values above
    # (and including) the matrix diagonal are incorporated into the newly
    # formed data frame. This piece of code drops the rows containing NA,
    # effectively retaining only those values *below* the matrix diagonal
    pwise_simil_df <- as.data.frame.table(
      x = pwise_simil_matrix,
      responseName = "value"
    ) %>%
      rename(
        population_1 = Var1,
        population_2 = Var2
      ) %>%
      drop_na() %>%
      # Remove numbers at the end of the population names
      mutate(
        population_1 = str_remove_all(
          string = population_1,
          pattern = "\\d"
        ),
        population_2 = str_remove_all(
          string = population_2,
          pattern = "\\d"
        )
      )


    ## First part of the summary output ##

    # Computes the mean pairwise similarity within the vancomycin group
    # and control group
    avg_similarity_df <- pwise_simil_df %>%
      group_by(population_1, population_2) %>%
      summarize(avg = mean(value)) %>%
      as_tibble() %>%
      filter(population_1 == population_2) %>%
      rename(
        group_1 = population_1,
        group_2 = population_2
      )


    ## Second part of the summary output ##

    difference <- avg_similarity_df[2, 3] - avg_similarity_df[1, 3]
    names(difference) <- "difference"

    # When the permutation argument is FALSE, then the calculate_similarity
    # function returns a list containing the *observed* mean pairwise
    # similarity within the vancomycin-adapted and control lines, and the
    # difference between these two values

    # When the permutation argument is TRUE, then the calculate_similarity
    # function returns a data frame containing the *permuted* mean pairwise
    # similarity within the vancomycin-adapted and control lines, and the
    # difference between these two values
    if (permutation == FALSE) {
      results <- list(
        avg_similarity_df,
        difference
      )
    } else {
      # Makes a data frame with three columns:
      # 1. Permuted average pairwise similarity - control populations
      # 2. Permuted average pairwise similarity - treated populations
      # 3. Difference between these permuted means
      results <- avg_similarity_df %>%
        mutate(
          group = ifelse(
            test = group_1 == "C",
            yes = "control",
            no = "treated"
          )
        ) %>%
        dplyr::select(
          -group_1,
          -group_2
        ) %>%
        pivot_wider(
          names_from = group,
          values_from = avg
        ) %>%
        mutate(difference = treated - control)
    }

    return(results)
  }

  # Function to perform permutation test

  perform_permutation_test <- function(.x, num_perm) {
    # Detect number of cores and register parallel backend
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    # Permuted distributions of the average pairwise similarity in
    # each treatment. By default, the number of permutations is set to
    # 10,000, but this can be changed based on user preferences
    results <- foreach(
      i = seq_len(num_perm),
      .combine = rbind,
      .packages = c("tidyverse", "proxy"),
      .export = "calculate_similarity"
    ) %dopar% {
      calculate_similarity(
        .x = .x,
        permutation = TRUE
      )
    }

    stopCluster(cl)

    return(results)
  }

  # Calculate *observed* average pairwise similarity within the
  # vancomycin-adapted and control lines
  similarity <- calculate_similarity(.x = .x)

  # Perform permutation test
  permutation_results <- perform_permutation_test(
    .x = .x,
    num_perm = num_perm
  )

  # Calculate significance
  significance <- sum(
    permutation_results$difference >= similarity[[2]]$difference
  ) / num_perm

  return(
    list(
      similarity = similarity,
      permutation_results = permutation_results,
      significance = significance
    )
  )
}


## -----------------------------------------------------------------------------
## Functions that (i) perform a multivariate logistic regression, (ii) test for
## coefficient estimate convergence, (iii) and plot cumulative mean and standard
## error of model coefficients.
## -----------------------------------------------------------------------------

## Regression function ##

perform_regression <- function(
    .x,
    iterations = 250,
    convergence_test = FALSE) {
  regression_df <- .x %>%
    as.data.frame() %>%
    rownames_to_column(var = "treatment") %>%
    # Remove numbers at the end of the population names; and,
    # if the value in the treatment column is T, replace it with 1.
    # Otherwise, replace it with 0
    mutate(
      treatment = str_remove_all(
        string = treatment,
        pattern = "\\d"
      ),
      treatment = ifelse(
        test = treatment == "T",
        yes = 1,
        no = 0
      )
    ) %>%
    # Remove columns in which the sum of all values is less than 3
    select_if(~ sum(.) >= 3)

  # Multivariate logistic regression routine

  # Detect number of cores and register parallel backend
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)

  coef_results <- foreach(
    i = seq_len(iterations),
    .combine = "rbind",
    .packages = c("glmnet")
  ) %dopar% {
    cvfit <- cv.glmnet(
      x = as.matrix(regression_df[, -1]),
      y = as.factor(regression_df$treatment),
      family = binomial(link = "probit"),
      gamma = seq(0, 1, 0.1),
      nfolds = 3,
      maxit = 1e7,
      relax = TRUE
    )

    fit_results <- glmnet(
      x = as.matrix(regression_df[, -1]),
      y = as.factor(regression_df$treatment),
      family = binomial(link = "probit"),
      lambda = cvfit$relaxed$lambda.min,
      alpha = cvfit$relaxed$gamma.min
    )

    coef_results <- as.matrix(fit_results$beta)

    return(coef_results)
  }

  stopCluster(cl)

  coef_results <- as.data.frame(coef_results)

  coef_results <- coef_results %>%
    rownames_to_column(var = "gene") %>%
    mutate(gene = sub("\\.\\d+$", "", gene)) %>%
    rename(coef = s0)

  # To test for convergence in coefficient estimates (i.e., how many times
  # do we need to run the regression to obtain stable estimates?) We'll track
  # the cumulative mean and standard error of the coefficients for each gene
  # over the iterations
  if (convergence_test == TRUE) {
    n_genes <- length(unique(coef_results$gene))
    n_iter <- nrow(coef_results) / n_genes

    coef_results <- coef_results %>%
      mutate(
        iteration = rep(
          x = 1:n_iter,
          each = n_genes
        )
      ) %>%
      pivot_wider(
        names_from = gene,
        values_from = coef
      )

    return(coef_results)
  }

  # For each gene, if the coefficient is non-zero less than 20% of the
  # number of total iterations, then exclude it. NOTE: we use 400 iterations
  # (convergence), so this function is set to that value by default
  # but it can be changed to a different value based on user preferences.
  idx <- coef_results %>%
    group_by(gene) %>%
    summarize(count = sum(coef != 0)) %>%
    filter(count >= (0.2 * iterations)) %>%
    pull(gene)

  coef_OR_df <- coef_results %>%
    filter(gene %in% idx) %>%
    # Exponentiate the coefficients to obtain the odds ratio
    mutate(OR = exp(coef))

  ## Determine the 95% confidence interval around the mean coefficient and OR
  ## for each gene
  regression_summary <- coef_OR_df %>%
    group_by(gene) %>%
    summarize(
      coef_mean = mean(coef),
      coef_SE = sd(coef) / sqrt(n()),
      OR_mean = mean(OR),
      OR_SE = sd(OR) / sqrt(n())
    ) %>%
    mutate(
      coef_lower_CI = coef_mean - (1.96 * coef_SE),
      coef_upper_CI = coef_mean + (1.96 * coef_SE)
    )

  idx_updated <- regression_summary %>%
    filter(abs(coef_mean) >= 0.05) %>%
    pull(gene)

  # Select columns of the regression_summary and regression_df data frames based
  # on the idx_updated vector (i.e., these coefficients were associated with
  # treatment status in at least 80 of the 400 iterations AND the absolute mean
  # coefficient is greater than or equal to 0.05)

  regression_summary <- regression_summary %>%
    filter(gene %in% idx_updated)

  coef_matrix <- as.matrix(regression_df[-1][, idx_updated])
  rownames(coef_matrix) <- rownames(.x)

  # Exponentiate the endpoints of the 95% CI for each coefficient to obtain
  # the confidence interval for the odds ratio
  regression_summary <- regression_summary %>%
    mutate(
      OR_lower_CI = exp(coef_lower_CI),
      OR_upper_CI = exp(coef_upper_CI)
    )

  # Order the columns in the coef_matrix by OR mean and list the
  # treated lines before the controls
  OR_matrix <- coef_matrix[
    c(
      paste0("T", 1:18),
      paste0("C", 1:87)
    ),
    order(regression_summary$OR_mean,
      decreasing = TRUE
    )
  ]

  results <- list(
    regression_summary,
    OR_matrix
  )

  return(results)
}

## Function to examine cumulative mean and standard error of model coefficients

perform_convergence_test <- function(.x) {
  cumulative_mean <- function(x) {
    sapply(
      X = 1:length(x),
      FUN = function(i) {
        mean(x[1:i])
      }
    )
  }

  cumulative_sd <- function(x) {
    sapply(
      X = 1:length(x),
      FUN = function(i) {
        if (i == 1) {
          return(NA)
        } else {
          sd(x[1:i])
        }
      }
    )
  }

  full_coef_df <- .x %>%
    pivot_longer(
      cols = -iteration,
      names_to = "gene",
      values_to = "coef"
    )

  # Compute cumulative statistics for each gene
  cumulative_stats <- full_coef_df %>%
    group_by(gene) %>%
    arrange(iteration) %>%
    mutate(
      cum_mean = cumulative_mean(coef),
      cum_sd   = cumulative_sd(coef),
      cum_n    = row_number(),
      cum_SE   = cum_sd / sqrt(cum_n)
    ) %>%
    ungroup()

  return(cumulative_stats)
}


## Generate convergence plots

generate_convergence_plots <- function(.x, .y) {
  # NOTE: We highlight those genes that have non-zero coefficients in at least
  # 20% of the iterations and have an absolute mean odds ratio (OR) >= 1.05.
  convergence_stats <- .x %>%
    mutate(
      highlight = ifelse(
        test = gene %in% .y[[1]]$gene,
        yes = "yes",
        no = "no"
      )
    )

  labels_df <- convergence_stats %>%
    filter(highlight == "yes") %>%
    group_by(gene) %>%
    filter(iteration == max(iteration)) %>%
    ungroup() %>%
    mutate(label_expr = sprintf("italic('%s')", gene))

  coef_cum_mean_plot <- convergence_stats %>%
    ggplot() +
    geom_line(
      data = convergence_stats %>%
        filter(highlight == "no"),
      aes(
        x = iteration,
        y = cum_mean,
        group = gene
      ),
      color = "grey80",
      linewidth = 0.5,
      alpha = 0.2
    ) +
    geom_line(
      data = convergence_stats %>%
        filter(highlight == "yes"),
      aes(
        x = iteration,
        y = cum_mean,
        group = gene
      ),
      color = "#606B91",
      linewidth = 1
    ) +
    ggrepel::geom_label_repel(
      data = labels_df,
      aes(
        x = iteration,
        y = cum_mean,
        label = label_expr
      ),
      parse = TRUE,
      nudge_x = 10,
      direction = "y",
      color = "#FFFFFF",
      fill = "#606B91",
      segment.color = "grey80",
      max.overlaps = Inf,
      force = 10,
      size = 2,
      box.padding = 0.1,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = 250,
      linetype = "dashed",
      color = "grey90"
    ) +
    scale_y_continuous(
      limits = c(-0.2, 1.2),
      breaks = seq(-0.2, 1.2, 0.2)
    ) +
    labs(
      x = "Iteration",
      y = "Cumulative Mean"
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(
        l = 13,
        r = 10,
        b = 25,
        t = 10,
        unit = "pt"
      )
    )

  coef_cum_se_plot <- convergence_stats %>%
    ggplot() +
    geom_line(
      data = convergence_stats %>%
        filter(highlight == "no"),
      aes(
        x = iteration,
        y = cum_SE,
        group = gene
      ),
      color = "grey80",
      linewidth = 0.5,
      alpha = 0.2
    ) +
    geom_line(
      data = convergence_stats %>%
        filter(highlight == "yes"),
      aes(
        x = iteration,
        y = cum_SE,
        group = gene
      ),
      color = "#606B91",
      linewidth = 1
    ) +
    geom_label_repel(
      data = labels_df,
      aes(
        x = iteration,
        y = cum_SE,
        label = label_expr
      ),
      parse = TRUE,
      nudge_x = 10,
      direction = "y",
      color = "#FFFFFF",
      fill = "#606B91",
      segment.color = "grey80",
      max.overlaps = Inf,
      force = 10,
      size = 2,
      box.padding = 0.1,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = 250,
      linetype = "dashed",
      color = "grey90"
    ) +
    scale_y_continuous(
      limits = c(0, 0.15)
    ) +
    labs(
      x = "Iterations",
      y = "Cumulative Standard Error"
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white")
    )

  convergence_plot <- plot_grid(
    coef_cum_mean_plot,
    coef_cum_se_plot,
    nrow = 2,
    labels = c("A", "B")
  )

  ggsave(
    filename = "convergence_plot.tiff",
    plot = convergence_plot,
    path = "Figures",
    device = "tiff",
    width = 6,
    height = 10,
    units = "in"
  )

  return(convergence_plot)
}

## Wrapper function
perform_regression_and_plots <- function(.x) {
  full_coef_df <- perform_regression(
    .x = .x,
    iterations = 1000,
    convergence_test = TRUE
  )

  convergence_stats <- perform_convergence_test(.x = full_coef_df)

  # The logistic regression model converges after ~250 iterations. We now re-run
  # it with 250 iterations to compute the odds ratios (ORs) for each gene.
  regression_results <- perform_regression(
    .x = binary_matrix,
    iterations = 250,
    convergence_test = FALSE
  )

  convergence_plots <- generate_convergence_plots(
    .x = convergence_stats,
    .y = regression_results
  )

  return(
    list(
      ... = regression_results,
      ... = convergence_plots
    )
  )
}


## -----------------------------------------------------------------------------
## Functions that plot odds ratios and a mutation heat map (Figure 2)
## -----------------------------------------------------------------------------

generate_OR_plot <- function(.x, .y) {
  OR_status <- .x %>%
    mutate(
      Group = ifelse(
        OR_mean > 1,
        "greater",
        "lesser"
      )
    )

  # Create the waterfall plot
  plot <- OR_status %>%
    ggplot(
      aes(
        x = factor(
          gene,
          levels = colnames(.y)
        ),
        y = OR_mean,
        color = Group
      )
    ) +
    geom_point(size = 3) +
    # geom_bar(
    #   stat = "identity",
    #   width = 0.9
    # ) +
    geom_errorbar(
      aes(
        ymin = OR_lower_CI,
        ymax = OR_upper_CI
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
    scale_y_continuous(limits = c(0.85, 3.25), breaks = seq(0.85, 3.25, 0.15)) +
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

generate_heatmap_plot <- function(.x, .y) {
  ## Initial data wrangling to prepare the mutation heat map

  # Identify genes with mean odds ratios greater than 1
  genes_OR_greater <- .x %>%
    filter(OR_mean > 1) %>%
    pull(gene)

  # Prepare data for heatmap
  OR_df <- .y %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Population",
      values_to = "Mutation"
    ) %>%
    mutate(
      OR_group = ifelse(
        test = Gene %in% genes_OR_greater,
        yes = "greater",
        no = "lesser"
      ),
      Treatment_group = ifelse(
        test = grepl(
          pattern = "T",
          x = Population
        ),
        yes = "Treated",
        no = "Control"
      ),
      Population = str_remove_all(
        string = Population,
        pattern = "\\D"
      )
    )

  # Heatmap plot
  heatmap_plot <- OR_df %>%
    ggplot() +
    facet_grid(
      factor(
        x = Treatment_group,
        levels = c("Treated", "Control")
      ) ~ .,
      scales = "free_y",
      space = "free_y",
      labeller = as_labeller(
        c(
          "Control" = "Control lines",
          "Treated" = "Vancomycin-adapted lines"
        )
      )
    ) +
    geom_tile(
      aes(
        y = fct_rev(
          factor(
            x = Population,
            levels = 1:87
          )
        ),
        x = factor(
          x = Gene,
          levels = colnames(.y)
        ),
        fill = Mutation
      ),
      data = filter(
        .data = OR_df,
        OR_group == "lesser"
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
            x = Population,
            levels = 1:87
          )
        ),
        x = factor(
          x = Gene,
          levels = colnames(.y)
        ),
        fill = Mutation
      ),
      data = filter(
        .data = OR_df,
        OR_group == "greater"
      ),
      color = "#EBEBEB"
    ) +
    scale_x_discrete(limits = colnames(.y)) +
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


## -----------------------------------------------------------------------------
## Function that performs Bayesian latent class analysis (BLCA) on the
## vancomycin-treated and control lines.
## -----------------------------------------------------------------------------

perform_blca <- function(.x, iterations = 1000, num_clusters = 5) {
  ## BLCA function ##

  blca <- function(OR_matrix, group_label) {
    browser()
    group_matrix <- OR_matrix[
      grep(
        pattern = group_label, # Group label corresponds to control or treated
        x = rownames(OR_matrix)
      ),
    ]

    # In the BLCA analysis, we only include genes that have a mutation in at
    # least 3 populations
    group_matrix <- group_matrix[, colSums(group_matrix) >= 3]

    # Subsample the group_matrix by randomly selecting 90% of the rows
    # and 90% of the columns
    local_matrix <- group_matrix[
      sample(
        seq_len(
          nrow(group_matrix)
        ),
        size = ceiling(
          nrow(group_matrix) * 0.9
        ),
        replace = FALSE
      ),
      sample(
        seq_len(
          ncol(group_matrix)
        ),
        size = ncol(group_matrix) * 0.9,
        replace = FALSE
      )
    ]

    tryCatch(
      expr = {
        # Main routine of LCA: given number of clusters (1:5) estimates the
        # model parameters and selects the best model based on AIC.
        local_blca_res <- map(
          .x = seq_len(num_clusters),
          ~ blca.em(
            X = local_matrix,
            G = .x,
            delta = 1,
            alpha = 0.5,
            beta = 0.5,
            restarts = 10,
            iter = 1000,
            start.vals = "across"
          )
        )

        best_k <- which.max(
          map_dbl(
            .x = local_blca_res,
            ~ .x$AIC
          )
        )

        best_res <- local_blca_res[[best_k]]

        # Cluster assignments are made and kept track of.
        if (best_k == 1) {
          cluster_assign <- setNames(
            rep(
              x = 1,
              times = nrow(local_matrix)
            ),
            rownames(local_matrix)
          )

          return(cluster_assign)
        } else {
          if ((nrow(best_res$Z) == nrow(local_matrix)) &
            is.null(
              rownames(best_res$Z)
            )
          ) {
            cluster_assign <- setNames(
              apply(
                X = best_res$Z,
                MARGIN = 1, # Function applied over rows
                FUN = which.max
              ),
              rownames(local_matrix)
            )

            return(cluster_assign)
          } else if (
            !is.null(
              rownames(best_res$Z)
            )
          ) {
            cluster_assign <- setNames(
              apply(
                X = best_res$Z,
                MARGIN = 1,
                FUN = which.max
              )
              [
                match(
                  apply(
                    X = local_matrix,
                    MARGIN = 1,
                    FUN = function(s) {
                      paste(s, collapse = "")
                    }
                  ),
                  table = rownames(best_res$Z)
                )
              ],
              rownames(local_matrix)
            )

            return(cluster_assign)
          } else {
            return(NULL)
          }
        }
      },
      error = function(cond) {
        return(NULL)
      },
      warning = function(cond) {
        return(NULL)
      }
    )
  }

  ## Parallelize the BLCA function ##

  parallelize_blca <- function(OR_matrix) {
    # Detect number of cores and register parallel backend
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    blca_res_em <- foreach(
      i = seq_len(iterations),
      .packages = c(
        "BayesLCA",
        "purrr"
      ),
      .export = "blca"
    ) %dopar% {
      map(
        # Factor over both the control and treated group labels
        .x = c(
          "^C",
          "^T"
        ),
        ~ blca(
          OR_matrix = OR_matrix,
          group_label = .x
        )
      )
    }

    stopCluster(cl)

    return(blca_res_em)
  }


  ## Function to construct the consensus matrices ##

  construct_frequency_matrix <- function(OR_matrix, blca_res_em, group_label) {
    group_matrix <- OR_matrix[
      grep(
        pattern = group_label, # Group label corresponds to control or treated
        x = rownames(OR_matrix)
      ),
    ]

    # In the BLCA analysis, we only include genes that have a mutation in at
    # least 3 populations
    group_matrix <- group_matrix[, colSums(group_matrix) >= 3]

    blca_res_em <- blca_res_em %>%
      map(
        ~ .x %>%
          discard(is.null)
      )

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
      message(
        sprintf(
          fmt = "Processing Run Index: %d",
          ... = i
        )
      )

      local_res <- blca_res_em[[i]]
      local_k <- max(local_res)

      local_clust <- as.matrix(
        Matrix::Diagonal(
          n = local_k + 1
        )
        [local_res, ]
      )

      local_clust <- local_clust %*% t(local_clust)

      # Reorder
      idx <- match(
        sorted_ids,
        table = names(local_res)
      )

      local_clust <- local_clust[idx, idx]
      local_clust[is.na(local_clust)] <- 0
      pair_matrix <- pair_matrix + local_clust

      local_clust <- as.matrix(
        Matrix::Diagonal(
          n = 2
        )
        [rep(
            x = 1,
            length = length(local_res)
          ), ]
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

    freq_matrix <- pair_matrix / (count_matrix + 1)
    freq_matrix <- Matrix::forceSymmetric(freq_matrix)
    freq_matrix <- as.matrix(freq_matrix)

    return(freq_matrix)
  }


  ## Hierarchical clustering and plotting function ##

  plot_clusters <- function(OR_matrix, freq_matrix, group_label) {
    group_matrix <- OR_matrix[
      grep(
        pattern = group_label, # Group label corresponds to control or treated
        x = rownames(OR_matrix)
      ),
    ]

    # In the BLCA analysis, we only include genes that have a mutation in at
    # least 3 populations
    group_matrix <- group_matrix[, colSums(group_matrix) >= 3]

    # Perform hierarchical clustering on the distance matrix
    dist_matrix <- as.dist(1.0 - freq_matrix)

    hc <- hclust(
      d = dist_matrix,
      method = "ward.D2"
    )

    cluster_colors <- setNames(
      RColorBrewer::brewer.pal(
        n = 7,
        name = "Paired"
      ),
      paste0("Cluster-", 1:7)
    )

    for (k in c(2, 3, 4, 5, 6, 7)) {
      clust_assign <- cutree(
        tree = hc,
        k = k
      )

      clust_res <- data.frame(
        "ID" = names(clust_assign),
        "Cluster" = paste0("Cluster-", clust_assign)
      )

      if (
        identical(
          x = freq_matrix,
          y = group_freq_matrices[[1]]
        )
      ) {
        output_tag <- "AIC_controlOnly"
      } else {
        output_tag <- "AIC_treatedOnly"
      }

      write.csv(
        x = clust_res,
        file = sprintf(
          fmt = "data\\cluster_results\\ClusterResults_%s_K%s.csv",
          ... = output_tag,
          ... = k
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
        breaks = seq(0, 1, 0.01),
        colors = colorRampPalette(
          colors = c("white", "grey30")
        )(
          length(
            seq(0, 1, 0.01)
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
        file = sprintf(
          fmt = "figures\\cluster_plots\\CC_Genomic_%s_K%s.pdf",
          ... = output_tag,
          ... = k
        ),
        width = 10,
        height = 6
      )

      draw(h_comb)
      dev.off()
    }
  }

  ## Function calls ##

  OR_matrix <- .x[[1]][[2]]

  blca_res_em <- parallelize_blca(OR_matrix = OR_matrix)

  # Separate the cluster assignments for the control and treated lines
  control_clusters <- blca_res_em %>%
    map(
      ~ .x[[1]]
    )

  treated_clusters <- blca_res_em %>%
    map(
      ~ .x[[2]]
    )

  # Construct the frequency matrix for the control and treated lines
  group_freq_matrices <- map2(
    .x = list(
      control_clusters,
      treated_clusters
    ),
    .y = list(
      "^C",
      "^T"
    ),
    ~ construct_frequency_matrix(
      OR_matrix = OR_matrix,
      blca_res_em = .x,
      group_label = .y
    )
  )

  cluster_plots <- map2(
    .x = list(
      group_freq_matrices[[1]],
      group_freq_matrices[[2]]
    ),
    .y = list(
      "^C",
      "^T"
    ),
    ~ plot_clusters(
      OR_matrix = OR_matrix,
      freq_matrix = .x,
      group_label = .y
    )
  )

  return(NULL)
}


## -----------------------------------------------------------------------------
## Function that calculates collateral response values.
## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
## Function that estimates summary stats of collateral response values.
## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
## Function to plot collateral responses (Figure 3C)
## -----------------------------------------------------------------------------

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

## -----------------------------------------------------------------------------
## Function that performs Mann-Whitney U tests on collateral response values.
## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
## Function that performs biserial correlation analysis on collateral response
## values and regulon membership (Figure 3)
## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
## Function that computes bootstrap collateral response scores (CRS) and
## percentile-based confidence intervals for each drug and mutational pathway
## and plots these distributions.
## -----------------------------------------------------------------------------

perform_bootstrap_and_plot_css <- function(.x, groups, iterations) {
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
        antibiotic_numeric = as.numeric(factor(antibiotic, levels = unique(antibiotic))),
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
        #alpha = 0.5,
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
          y = antibiotic,         # use factor for centering
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
      # #scale_y_discrete(
      #   breaks = unique(boot_res_both$antibiotic_numeric),
      #   labels = unique(boot_res_both$antibiotic),
      #   name = "Antibiotic"
      # ) +
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