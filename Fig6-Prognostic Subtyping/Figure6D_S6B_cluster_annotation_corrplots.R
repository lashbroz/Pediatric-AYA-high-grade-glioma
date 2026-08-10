#!/usr/bin/env Rscript

# ==============================================================================
# File: Figure6D_S6B_cluster_annotation_corrplots.R
# ==============================================================================
#
# Title:
#   Figure 6D and Figure S6B —
#   enrichment of molecular, clinical, and reference subtype annotations
#   across proteomic subtype clusters
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   This script generates corrplot-style enrichment visualizations evaluating
#   associations between collapsed proteomic subtype clusters and reference
#   subtype, clinical, and molecular annotations in the Discovery cohort.
#
#   Sex-stratified proteomic subtypes are collapsed as:
#
#       M1/F1 -> Cluster 1
#       M2/F2 -> Cluster 2
#       M3/F3 -> Cluster 3
#
#   Figure 6D:
#     - H3-3A mutation enrichment across proteomic subtype clusters
#
#   Figure S6B:
#     - Age class enrichment
#     - Cancer group enrichment
#     - Wang classification enrichment
#     - Garofano classification enrichment
#
# Statistical framework:
#   - One-sided Fisher exact tests were performed independently for each
#     annotation-cluster combination to evaluate enrichment relative to the
#     remainder of the contingency table.
#   - P-values were adjusted across all comparisons within each contingency
#     matrix using the Benjamini-Hochberg method.
#   - Significant enrichments are annotated directly on corrplots using
#     FDR-adjusted significance thresholds.
#
# Input files:
#   - data/STable1.xlsx
#   - data/STable6.xlsx
#   - data/cDisc_mutation_10192023.tsv
#
# Output files:
#
#   Figure 6D:
#     - Figure6D_cluster_annotation_corrplot_H3_3A_mutation_status.pdf
#
#   Figure S6B:
#     - FigureS6B_cluster_annotation_corrplot_Age_class.pdf
#     - FigureS6B_cluster_annotation_corrplot_Cancer_group.pdf
#     - FigureS6B_cluster_annotation_corrplot_Wang_class.pdf
#     - FigureS6B_cluster_annotation_corrplot_Garofano_class.pdf
#
# ==============================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(corrplot)
  library(viridis)
})

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

clinical_file <- "../data/STable1.xlsx"
subtype_file  <- "../data/STable6.xlsx"
mutation_file <- "../data/cDisc_mutation_10192023.tsv"

# ------------------------------------------------------------------------------
# Read input data
# ------------------------------------------------------------------------------

clinical_df <- read_excel(clinical_file, sheet = "ClinicalTable") %>%
  as.data.frame()

subtype_df <- read_excel(subtype_file, sheet = "Subtype-cDisc") %>%
  as.data.frame()

mutation_df <- read.delim(
  mutation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

clean_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "", "NaN")] <- NA
  x
}

extract_h3_3a_status <- function(mutation_df, sample_ids) {
  h3_row <- mutation_df[mutation_df$ApprovedGeneSymbol == "H3-3A", ]
  
  if (nrow(h3_row) != 1) {
    stop("Expected exactly one H3-3A row; found ", nrow(h3_row))
  }
  
  matched_cols <- match(sample_ids, colnames(mutation_df))
  
  h3_values <- rep(NA_real_, length(sample_ids))
  matched <- !is.na(matched_cols)
  
  h3_values[matched] <- as.numeric(
    as.matrix(h3_row[, matched_cols[matched], drop = FALSE])
  )
  
  h3_values[is.na(h3_values)] <- 0
  h3_values
}

prepare_analysis_data <- function(clinical_df, subtype_df, mutation_df) {
  data <- clinical_df
  
  data$id <- as.character(data$id)
  subtype_df$id <- as.character(subtype_df$id)
  
  subtype_match <- match(data$id, subtype_df$id)
  
  data$wang.cl <- clean_na(subtype_df$wang.cl[subtype_match])
  data$garofano.cl <- clean_na(subtype_df$garofano.cl[subtype_match])
  
  data$cl <- subtype_df$protein.subtype[subtype_match]
  data$cl <- gsub("(M|F)", "", data$cl)
  data$cl <- factor(data$cl, levels = c("1", "2", "3"))
  
  data$H3.3A <- extract_h3_3a_status(
    mutation_df = mutation_df,
    sample_ids = data$id
  )
  
  data$H3.3A <- factor(
    data$H3.3A,
    levels = c(0, 1),
    labels = c("WT", "Mut")
  )
  
  data$age.class <- factor(
    data$cDisc_age_class_name_derived,
    levels = c("PED", "ADO", "YA", "ADULT")
  )
  
  data$cancer.group <- clean_na(data$Disc_Cancer_Group)
  
  data %>%
    filter(!is.na(cl))
}

make_contingency_table <- function(data, var) {
  plot_data <- data %>%
    filter(
      !is.na(cl),
      !is.na(.data[[var]])
    )
  
  tab <- table(
    plot_data[[var]],
    plot_data$cl
  )
  
  tab <- as.matrix(tab)
  colnames(tab) <- paste0("C", colnames(tab))
  
  tab
}

run_fisher_enrichment <- function(tab, adjust_method = "BH") {
  p_mat <- matrix(
    1,
    nrow = nrow(tab),
    ncol = ncol(tab),
    dimnames = dimnames(tab)
  )
  
  odds_mat <- matrix(
    NA_real_,
    nrow = nrow(tab),
    ncol = ncol(tab),
    dimnames = dimnames(tab)
  )
  
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      a <- tab[i, j]
      b <- sum(tab[i, ]) - a
      c <- sum(tab[, j]) - a
      d <- sum(tab) - a - b - c
      
      fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
      
      fit <- fisher.test(
        fisher_mat,
        alternative = "greater"
      )
      
      p_mat[i, j] <- fit$p.value
      odds_mat[i, j] <- unname(fit$estimate)
    }
  }
  
  fdr_mat <- matrix(
    p.adjust(as.vector(p_mat), method = adjust_method),
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )
  
  list(
    p_value = p_mat,
    fdr = fdr_mat,
    odds_ratio = odds_mat
  )
}

make_stars <- function(fdr_mat) {
  stars <- matrix(
    "",
    nrow = nrow(fdr_mat),
    ncol = ncol(fdr_mat),
    dimnames = dimnames(fdr_mat)
  )
  
  stars[fdr_mat < 0.05] <- "*"
  stars[fdr_mat < 0.01] <- "**"
  stars[fdr_mat < 0.001] <- "***"
  
  stars
}

save_corrplot <- function(tab,
                          fdr_mat,
                          output_pdf,
                          plot_title = NULL,
                          width = 8,
                          height = 6) {
  stars <- make_stars(fdr_mat)
  
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = TRUE)) else getwd()
  output_dir <- file.path(script_dir, "output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pdf(file.path(output_dir, output_pdf), width = width, height = height, useDingbats = FALSE)
  
 corrplot(
  tab,
  is.corr = FALSE,
  method = "circle",
  col = colorRampPalette(
    viridis(3, option = "viridis", alpha = 0.75, direction = -1)
  )(75),
  cl.ratio = .5,
  cl.offset = 3,
  tl.col = "black",
  mar = c(1, 2, 3, 1)
)
  
  text(
    x = rep(seq_len(ncol(tab)), each = nrow(tab)),
    y = rep(nrow(tab):1, times = ncol(tab)),
    labels = stars,
    col = "white",
    cex = 1.6,
    font = 2
  )
  
  dev.off()
  
  message("Wrote: ", output_pdf)
}

run_annotation_corrplot <- function(data,
                                    var,
                                    label,
                                    figure_prefix,
                                    adjust_method = "BH") {
  message("Running annotation enrichment for: ", label)
  
  tab <- make_contingency_table(
    data = data,
    var = var
  )
  
  stats <- run_fisher_enrichment(
    tab = tab,
    adjust_method = adjust_method
  )
  
  safe_label <- gsub("[^A-Za-z0-9]+", "_", label)
  
  output_prefix <- paste0(
    figure_prefix,
    "_cluster_annotation_corrplot_",
    safe_label
  )
  
  save_corrplot(
    tab = tab,
    fdr_mat = stats$fdr,
    output_pdf = paste0(output_prefix, ".pdf"),
    #plot_title = label,
    width = 8,
    height = 6
  )
  
  invisible(
    list(
      variable = var,
      label = label,
      counts = tab,
      p_value = stats$p_value,
      fdr = stats$fdr,
      odds_ratio = stats$odds_ratio
    )
  )
}

# ------------------------------------------------------------------------------
# Prepare merged analysis data
# ------------------------------------------------------------------------------

data <- prepare_analysis_data(
  clinical_df = clinical_df,
  subtype_df = subtype_df,
  mutation_df = mutation_df
)

# ------------------------------------------------------------------------------
# Variables to analyze
# ------------------------------------------------------------------------------

variables_to_plot <- list(
  list(
    var = "H3.3A",
    label = "H3_3A_mutation_status",
    figure_prefix = "Figure6D"
  ),
  list(
    var = "age.class",
    label = "Age_class",
    figure_prefix = "FigureS6B"
  ),
  list(
    var = "cancer.group",
    label = "Cancer_group",
    figure_prefix = "FigureS6B"
  ),
  list(
    var = "wang.cl",
    label = "Wang_class",
    figure_prefix = "FigureS6B"
  ),
  list(
    var = "garofano.cl",
    label = "Garofano_class",
    figure_prefix = "FigureS6B"
  )
)

# ------------------------------------------------------------------------------
# Run analyses
# ------------------------------------------------------------------------------

results <- list()

for (item in variables_to_plot) {
  results[[item$label]] <- run_annotation_corrplot(
    data = data,
    var = item$var,
    label = item$label,
    figure_prefix = item$figure_prefix,
    adjust_method = "BH"
  )
}

message("Done.")
