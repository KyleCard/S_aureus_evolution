### Functions associated with the Card and Crozier et al. (2025) R Notebook
### Authors: Kyle Card and Arda Durmaz


## -----------------------------------------------------------------------------
## Function that plots vancomycin MIC data and generates Figure 1 (with
## the experimental schema).
## -----------------------------------------------------------------------------

generate_vanc_plot <- function(.x) {
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
  # mutations object and populates the matrix with 0s.
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
## Function that performs similarity and bootstrap analyses on the qualifying
## gene-level mutations in the vancomycin-treated and control lines.
## -----------------------------------------------------------------------------

perform_genomic_analysis <- function(.x, iterations = 10000) {
  calculate_similarity <- function(.x, bootstrap = FALSE) {
    if (bootstrap == TRUE) {
      # Sample the lines, with replacement, 105 times.
      .x <- .x[
        sample(
          nrow(.x),
          replace = TRUE
        ),
      ]

      # Rename the rows of the sampled dataset - some rows are sampled
      # more than once. By default R will append a number to the end of the
      # row name to ensure unique names. These rows are then treated separately
      # in the downstream analyses. This code chunk renames the rows to
      # C1, C2, ..., C87 and T1, T2, ..., T18 etc to avoid this issue.
      rownames(.x) <- c(
        paste0("C", 1:87),
        paste0("T", 1:18)
      )
    }

    # Computes Dice's Similarity Coefficient for all possible replicate pairs.
    # Populates a matrix with these values.
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
    # NA to ease downstream wrangling.
    diag(pwise_simil_matrix) <- NA
    pwise_simil_matrix[
      upper.tri(pwise_simil_matrix)
    ] <- NA

    # Converted matrix to a data frame for analysis. NA values above
    # (and including) the matrix diagonal are incorporated into the newly
    # formed data frame. This piece of code drops the rows containing NA,
    # effectively retaining only those values *below* the matrix diagonal.
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
    # and control group.
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

    # When the bootstrap argument is FALSE, then the calculate_similarity
    # function returns a list containing the *observed* mean pairwise
    # similarity within the vancomycin-adapted and control lines, and the
    # difference between these two values.

    # When the bootstrap argument is TRUE, then the calculate_similarity
    # function returns a data frame containing the *bootstrap* mean pairwise
    # similarity within the vancomycin-adapted and control lines, and the
    # difference between these two values.
    if (bootstrap == FALSE) {
      results <- list(
        avg_similarity_df,
        difference
      )
    } else {
      # Makes a data frame with three columns:
      # 1. Bootstrap average pairwise similarity - control populations
      # 2. Bootstrap average pairwise similarity - treated populations
      # 3. Difference between these bootstrap means
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

  # Function to perform bootstrap analysis
  perform_bootstrap_analysis <- function(.x, iterations) {
    # Detect number of cores and register parallel backend
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    # Bootstrap distributions of the average pairwise similarity in
    # each treatment
    results <- foreach(
      i = seq_len(iterations),
      .combine = rbind,
      .packages = c("tidyverse", "proxy"),
      .export = "calculate_similarity"
    ) %dopar% {
      calculate_similarity(
        .x = .x,
        bootstrap = TRUE
      )
    }

    stopCluster(cl)

    return(results)
  }

  # Calculate *observed* average pairwise similarity within the
  # vancomycin-adapted and control lines
  similarity <- calculate_similarity(.x = .x)

  # Perform bootstrap analysis
  bootstrap_results <- perform_bootstrap_analysis(
    .x = .x,
    iterations = iterations
  )

  # Calculate significance
  significance <- sum(
    bootstrap_results$difference >= similarity[[2]]$difference
  ) / iterations

  return(
    list(
      similarity = similarity,
      bootstrap_results = bootstrap_results,
      significance = significance
    )
  )
}


## -----------------------------------------------------------------------------
## Function that performs a multivariate logistic regression and Bayesian
## latent class analysis (BLCA) on the vancomycin-treated and control lines.
## -----------------------------------------------------------------------------

perform_regression_blca <- function(.x, iterations = 1000, num_clusters = 5) {
  ## Regression function ##

  perform_regression <- function(.x) {
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

    # Multivariate logistic regression routine over 1,000 iterations to select
    # robust non-zero coefficients

    # Detect number of cores and register parallel backend
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    coef_res <- foreach(
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

      fit_res <- glmnet(
        x = as.matrix(regression_df[, -1]),
        y = as.factor(regression_df$treatment),
        family = binomial(link = "probit"),
        lambda = cvfit$relaxed$lambda.min,
        alpha = cvfit$relaxed$gamma.min
      )

      coef_res <- as.matrix(fit_res$beta)

      return(coef_res)
    }

    stopCluster(cl)

    coef_res <- as.data.frame(coef_res)

    coef_res <- coef_res %>%
      rownames_to_column(var = "genes") %>%
      rename(coef = s0)

    # Remove the decimal and number from the end of the gene names
    coef_res$genes <- str_remove_all(
      string = coef_res$genes,
      pattern = "\\..*"
    )

    # For each gene, if the coefficient is non-zero less than 20% of the
    # number of total iterations, then exclude it. NOTE: in the manuscript, we
    # used 1,000 iterations, so this function is set to that value by default
    # but it can be changed to a different value based on user preferences.
    idx <- coef_res %>%
      group_by(genes) %>%
      summarize(count = sum(coef != 0)) %>%
      filter(count >= (0.2 * iterations)) %>%
      pull(genes)

    coef_OR_filtered <- coef_res %>%
      filter(genes %in% idx) %>%
      # Exponentiate the coefficients to obtain the odds ratio
      mutate(OR = exp(coef))

    ## Determine the 95% confidence interval around the mean coefficient and OR
    ## for each gene
    regression_summary <- coef_OR_filtered %>%
      group_by(genes) %>%
      summarize(
        coef_mean = mean(coef),
        coef_SE = sd(coef) / sqrt(n()),
        OR_mean = mean(OR),
        OR_SE = sd(OR) / sqrt(n())
      ) %>%
      mutate(
        coef_lower_CI = coef_mean - (1.96 * coef_SE),
        coef_upper_CI = coef_mean + (1.96 * coef_SE),
        OR_lower_CI = OR_mean - (1.96 * OR_SE),
        OR_upper_CI = OR_mean + (1.96 * OR_SE)
      )

    # Select columns of the regression_df data frame based on the idx vector
    # (i.e., these coefficients were associated with treatment status in at
    # least 200 of the 1,000 iterations)
    coef_matrix <- as.matrix(regression_df[-1][, idx])
    rownames(coef_matrix) <- rownames(.x)

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


  ## BLCA function ##

  perform_blca <- function(OR_matrix, group_label) {
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
      .export = "perform_blca"
    ) %dopar% {
      map(
        # Factor over both the control and treated group labels
        .x = c(
          "^C",
          "^T"
        ),
        ~ perform_blca(
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

  regression_results <- perform_regression(.x = binary_matrix)

  blca_res_em <- parallelize_blca(OR_matrix = regression_results[[2]])

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
      OR_matrix = regression_results[[2]],
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
      OR_matrix = regression_results[[2]],
      freq_matrix = .x,
      group_label = .y
    )
  )

  return(regression_results)
}


## -----------------------------------------------------------------------------
## Function that plots odds ratios and a mutation heat map.
## -----------------------------------------------------------------------------

generate_OR_heatmap_plot <- function(.x, .y) {
  OR_status <- .x %>%
    mutate(
      Group = ifelse(
        OR_mean > 1,
        "greater",
        "lesser"
      )
    )

  # Create the waterfall plot
  waterfall_plot <- OR_status %>%
    ggplot(
      aes(
        y = factor(
          genes,
          levels = rev(
            colnames(.y)
          )
        ),
        x = OR_mean,
        fill = Group
      )
    ) +
    geom_bar(
      stat = "identity",
      width = 0.9
    ) +
    geom_errorbarh(
      aes(
        xmin = OR_lower_CI,
        xmax = OR_upper_CI
      ),
      height = 0.2
    ) +
    ylab("Gene") +
    xlab("Average odds ratio") +
    scale_fill_manual(values = c("#606B91", "#E0E0E0")) +
    theme_cowplot() +
    theme(
      axis.text.y = element_text(
        angle = 0,
        face = "italic"
      ),
      panel.grid.major = element_line(color = "#EBEBEB"),
      plot.margin = margin(
        t = 25.5,
        b = 9
      ),
      legend.position = "none"
    )

  ## Initial data wrangling to prepare the mutation heat map

  # Genes with mean odds ratios greater than 1
  genes_OR_greater <- .x %>%
    filter(OR_mean > 1) %>%
    pull(genes)

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
        Gene %in% genes_OR_greater,
        yes = "greater",
        no = "lesser"
      ),
      Treatment_group = ifelse(
        grepl(
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

  heatmap_plot <- OR_df %>%
    ggplot() +
    facet_grid(
      ~ factor(
        Treatment_group,
        levels = c("Treated", "Control")
      ),
      scales = "free_x",
      space = "free_x",
      labeller = as_labeller(
        c(
          "Treated" = "Vancomycin-adapted lines",
          "Control" = "Control lines"
        )
      )
    ) +
    # Apply a grey (#E7E5DF) fill to the tiles with an odds ratio less than 1
    geom_tile(
      aes(
        x = factor(
          Population,
          levels = c(1:87)
        ),
        y = factor(
          Gene,
          levels = rev(
            colnames(.y)
          )
        ),
        fill = Mutation
      ),
      filter(
        .data = OR_df,
        OR_group == "lesser"
      ),
      color = "#EBEBEB"
    ) +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#E0E0E0"
    ) +

    # Start a new scale
    new_scale_fill() +

    # Apply a blue (#606B91) fill to the tiles with an odds ratio greater
    # than 1
    geom_tile(
      aes(
        x = factor(
          Population,
          levels = c(1:87)
        ),
        y = factor(
          Gene,
          levels = rev(
            colnames(.y)
          )
        ),
        fill = Mutation
      ),
      filter(
        .data = OR_df,
        OR_group == "greater"
      ),
      color = "#EBEBEB"
    ) +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#606B91"
    ) +
    xlab("Population") +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  combined_plot <- plot_grid(
    ... = waterfall_plot,
    ... = heatmap_plot,
    nrow = 1,
    ncol = 2,
    rel_widths = c(0.25, 1)
  ) +
    theme(plot.background = element_rect(fill = "white"))

  return(combined_plot)
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
## Function that calculates the probability of exhibiting a drug response.
## -----------------------------------------------------------------------------

calculate_likelihood <- function(.x) {
  .x <- .x %>%
    filter(antibiotic != "VAN")

  prob_vec <- c()

  # For each antibiotic in the summary_stats data frame, we calculate the
  # median collateral response values using the mean_CR column.
  median_df <- .x %>%
    group_by(antibiotic) %>%
    summarize(median = median(mean_CR))

  # We iterate through the unique antibiotics in the summary_stats
  # data frame and calculates the probability of a population exhibiting a drug
  # response. This probability is defined as the fraction of evolved populations
  # that have a collateral response value that has the same sign as the median
  # collateral response value. The probability of a population exhibiting a drug
  # response is calculated for each antibiotic and stored in the prob_vec vector
  for (i in unique(.x$antibiotic)) {
    L <- 0

    ab_frame <- .x %>%
      filter(antibiotic == i)

    med_CR <- ab_frame %>%
      pull(mean_CR) %>%
      median()

    for (j in ab_frame$mean_CR) {
      if (sign(j) == sign(med_CR)) {
        # We initialize a counter (L) at zero and increment it by a fraction
        # (1 / total number of populations) for each mean_CR value that is in
        # concordance with the median.
        L <- L + (1 / nrow(ab_frame))
      }
    }

    prob_vec <- c(prob_vec, L)
  }

  probability_df <- tibble(
    antibiotic = unique(.x$antibiotic),
    probability = prob_vec
  )

  # Join the median_df and probability_df objects by antibiotic
  df <- median_df %>%
    left_join(
      y = probability_df,
      by = "antibiotic"
    )

  return(df)
}


## -----------------------------------------------------------------------------
## Function to plot collateral responses and likelihoods
## -----------------------------------------------------------------------------

generate_CR_likelihood_plot <- function(.x, .y) {
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

  .x <- .x %>%
    filter(antibiotic != "VAN")

  CR_boxplot <- .x %>%
    ggplot(
      aes(
        x = antibiotic,
        y = mean_CR,
        fill = antibiotic
      )
    ) +
    geom_boxplot() +
    geom_jitter(
      shape = 16,
      position = position_jitter(0.2),
      size = 2
    ) +
    xlab("Antibiotic") +
    ylab(
      expression(
        "Log"[2] ~ "MIC"["evolved"] ~ "-" ~ "Log"[2] ~ "MIC"["ancestor"]
      )
    ) +
    scale_fill_manual(
      values = c(
        "#D09384", "#E4BDB2", "#F7E8E6", "#FFFFFF",
        "#F3EEF9", "#E3DAF3", "#D4C3EA", "#CAB7E6"
      )
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white")
    )

  likelihood_plot <- .y %>%
    ggplot(
      aes(
        x = probability,
        y = median
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_point(
      size = 5.5,
      fill = c(
        "#E4BDB2", "#F3EEF9", "#D09384", "#CAB7E6",
        "#E3DAF3", "#D4C3EA", "#FFFFFF"
      ),
      pch = 21
    ) +
    # Label each point with its corresponding antibiotic
    geom_text(
      aes(label = antibiotic),
      hjust = -0.5,
      vjust = -0.5,
      size = 5
    ) +
    scale_x_continuous(
      limits = c(0.5, 1),
      breaks = seq(0.5, 1, 0.1)
    ) +
    scale_y_continuous(
      limits = c(-1, 2),
      breaks = seq(-1, 2, 0.5)
    ) +
    xlab("Likelihood of drug response") +
    ylab(
      expression(
        "Log"[2] ~ "MIC"["evolved"] ~ "-" ~ "Log"[2] ~ "MIC"["ancestor"]
      )
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white")
    )

  combined_plot <- plot_grid(
    ... = CR_boxplot,
    ... = likelihood_plot,
    ncol = 2,
    rel_heights = c(0.5, 0.5),
    labels = c("A", "B"),
    label_size = 16
  )

  return(combined_plot)
}


## -----------------------------------------------------------------------------
## Function to plot association between cluster assignment and CR values
## -----------------------------------------------------------------------------

generate_CR_cluster_plot <- function(.x, clusters = 2) {
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

  cluster_assignments <- read_csv(
    file = paste0(
      ... = "data/cluster_results/ClusterResults_AIC_treatedOnly_K",
      ... = clusters,
      ... = ".csv"
    ),
    show_col_types = FALSE
  )

  # Assign clusters to the populations (needs perform_regression_blca to be
  # run first). The cluster assignment data is stored in the file
  # "ClusterResults_AIC_treatedOnly_K$.csv", where $ is the number of clusters
  # (K = 2 through 7).
  cluster_assignments <- cluster_assignments %>%
    mutate(
      population = str_remove(
        string = ID,
        pattern = "T"
      )
    ) %>%
    dplyr::select(-ID)

  CR_cluster_data <- .x %>%
    left_join(
      cluster_assignments,
      by = c("population" = "population")
    )

  CR_cluster_plot <- CR_cluster_data %>%
    ggplot(
      aes(
        x = antibiotic,
        y = mean_CR,
        fill = Cluster
      )
    ) +
    geom_boxplot() +
    xlab("Antibiotic") +
    ylab(
      expression(
        "Log"[2] ~ "MIC"["evolved"] ~ "-" ~ "Log"[2] ~ "MIC"["ancestor"]
      )
    ) +
    scale_fill_manual(
      values = c("#4ABCD7", "#E74B34")
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white")
    )

  return(CR_cluster_plot)
}
