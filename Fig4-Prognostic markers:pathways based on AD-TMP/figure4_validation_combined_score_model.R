#!/usr/bin/env Rscript

# ============================================================
# Figure 4 / Figure S4 | Validation of Male-Derived MORC2 and PHKG2 Prognostic Markers
# File: Figure4_validation_combined_score_model.R
#
# Description:
#   Performs validation survival analysis for MORC2 and PHKG2 in the
#   independent validation cohort.
#
#   MORC2 and PHKG2 were selected from the male-specific cDiscovery
#   prognostic model and are evaluated here as frozen male-derived
#   prognostic markers.
#
#   The script:
#     - Loads validation clinical and protein data from STable1.xlsx
#     - Prepares OS survival object
#     - Adds clinical adjustment covariates
#     - Tests MORC2 and PHKG2 as continuous markers
#     - Generates median-dichotomized adjusted survival plots
#     - Outputs HR tables and survival PDFs
#
# Inputs:
#   - ../data/STable1.xlsx
#
# Outputs:
#   Written to the current figure folder.
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================



suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(readxl)
})

set.seed(123)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

data_type <- "Validation"

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)))
  }
  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop("Cannot determine script directory for `", script_name, "`. Run from the script folder or with Rscript.", call. = FALSE)
}
script_dir <- resolve_script_dir("Figure4_validation_combined_score_model.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sex_colors <- c(
  Male = "#0707CF",
  Female = "#CC0303",
  Both = "black"
)

marker_genes <- c("MORC2", "PHKG2")

text_size <- 13

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

safe_scale_rows <- function(x) {
  x <- as.matrix(x)
  t(scale(t(x)))
}

make_adjustment_data <- function(adjustment_covariates, n_rows) {
  data.frame(
    matrix(
      0,
      nrow = n_rows,
      ncol = length(adjustment_covariates),
      dimnames = list(NULL, adjustment_covariates)
    )
  )
}

make_median_group <- function(x) {
  cut(
    x,
    breaks = c(
      min(x, na.rm = TRUE),
      median(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    ),
    include.lowest = TRUE,
    labels = c("Low", "High")
  )
}

save_pdf <- function(plot_object, filename, width, height) {
  pdf(file.path(output_dir, filename), width = width, height = height)
  print(plot_object)
  dev.off()
}

# ------------------------------------------------------------------------------
# Load validation data
# ------------------------------------------------------------------------------

vclinical <- readxl::read_xlsx(
  "../data/STable1.xlsx",
  sheet = 5,
  na = "NA"
)

vprotein <- readxl::read_xlsx(
  "../data/STable1.xlsx",
  sheet = 6,
  na = "NA"
)

vclinical <- as.data.frame(vclinical)
vprotein <- as.data.frame(vprotein)

# ------------------------------------------------------------------------------
# Prepare validation clinical data
# ------------------------------------------------------------------------------

vclinical <- vclinical %>%
  filter(
    age_class_derived %in% c("PED", "ADO"),
    Gender %in% c("Male", "Female"),
    surv_include
  ) %>%
  mutate(
    SurvObj = survival::Surv(survival_days, os_status),
    tumor = 1 - is_initial_tumor,
    ado = as.numeric(age_class_derived == "ADO"),
    Male = as.numeric(Gender == "Male"),
    mut_H33A = H33A_mut
  )

# Harmonize location/grade covariates if present under common names
if (!"Midline" %in% colnames(vclinical) && "tumor_loc" %in% colnames(vclinical)) {
  vclinical$Midline <- as.numeric(vclinical$tumor_loc == "Midline")
}

if (!"Grade" %in% colnames(vclinical) && "WHO_Grade" %in% colnames(vclinical)) {
  vclinical$Grade <- as.numeric(vclinical$WHO_Grade %in% 3)
}

# ------------------------------------------------------------------------------
# Prepare validation protein matrix
# ------------------------------------------------------------------------------

vprotein_matrix <- vprotein
rownames(vprotein_matrix) <- vprotein_matrix$Gene

vprotein_matrix <- vprotein_matrix[
  ,
  !colnames(vprotein_matrix) %in% "Gene",
  drop = FALSE
]

colnames(vprotein_matrix) <- gsub(
  "\\.",
  "-",
  gsub("^X", "", colnames(vprotein_matrix))
)

available_marker_genes <- intersect(marker_genes, rownames(vprotein_matrix))
missing_marker_genes <- setdiff(marker_genes, available_marker_genes)

if (length(missing_marker_genes) > 0) {
  message(
    "Marker genes missing from validation protein matrix: ",
    paste(missing_marker_genes, collapse = ", ")
  )
}

if (length(available_marker_genes) == 0) {
  stop("None of MORC2/PHKG2 were found in the validation protein matrix.")
}

validation_ids <- intersect(vclinical$id, colnames(vprotein_matrix))

vclinical <- vclinical[
  vclinical$id %in% validation_ids,
  ,
  drop = FALSE
]

vprotein_sub <- vprotein_matrix[
  available_marker_genes,
  vclinical$id,
  drop = FALSE
]

vprotein_sub <- safe_scale_rows(vprotein_sub)
vprotein_sub <- t(vprotein_sub)

validation_data <- data.frame(vclinical, vprotein_sub)

# ------------------------------------------------------------------------------
# Define adjustment covariates
# ------------------------------------------------------------------------------

adjustment_covariates <- c(
  "mut_H33A",
  "ado",
  "Midline",
  "Grade",
  "tumor",
  "Male"
)

adjustment_covariates <- intersect(
  adjustment_covariates,
  colnames(validation_data)
)

# ------------------------------------------------------------------------------
# Single-marker validation function
# ------------------------------------------------------------------------------

plot_single_marker_validation <- function(pred, data, adjustment_covariates) {
  
  if (!pred %in% colnames(data)) {
    message("Skipping ", pred, ": not found in validation_data")
    return(NULL)
  }
  
  model_covariates <- c(pred, adjustment_covariates)
  
  data <- data[
    complete.cases(data[, model_covariates, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  continuous_formula <- as.formula(
    paste0(
      "SurvObj ~ ",
      pred,
      if (length(adjustment_covariates) > 0) {
        paste0(" + ", paste(adjustment_covariates, collapse = " + "))
      } else {
        ""
      }
    )
  )
  
  null_formula <- as.formula(
    paste0(
      "SurvObj ~ ",
      paste(adjustment_covariates, collapse = " + ")
    )
  )
  
  binary_formula <- as.formula(
    paste0(
      "SurvObj ~ cl",
      if (length(adjustment_covariates) > 0) {
        paste0(" + ", paste(adjustment_covariates, collapse = " + "))
      } else {
        ""
      }
    )
  )
  
  fit_continuous <- coxph(
    continuous_formula,
    control = coxph.control(iter.max = 2000),
    data = data
  )
  
  fit_null <- coxph(
    null_formula,
    control = coxph.control(iter.max = 2000),
    data = data
  )
  
  anova_p <- anova(fit_continuous, fit_null)[2, 4]
  
  data$cl <- make_median_group(data[, pred])
  data$cl <- factor(data$cl, levels = c("Low", "High"))
  
  fit_binary <- coxph(
    binary_formula,
    control = coxph.control(iter.max = 2000),
    data = data
  )
  
  print(pred)
  print(summary(fit_continuous))
  print(summary(fit_binary))
  
  # --------------------------------------------------------------------------
  # HR table for continuous marker model
  # --------------------------------------------------------------------------
  
  marker_ci <- exp(confint(fit_continuous, level = 0.90))
  
  marker_hr_table <- data.frame(
    Variable = rownames(marker_ci),
    HR = exp(coef(fit_continuous)),
    CI_lower = marker_ci[, 1],
    CI_upper = marker_ci[, 2],
    stringsAsFactors = FALSE
  )
  
  write.table(
    marker_hr_table,
    file = file.path(
      output_dir,
      paste0(pred, "_MaleDerived_", data_type, "_hr_table.tsv")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
  
  # --------------------------------------------------------------------------
  # HR plot
  # --------------------------------------------------------------------------
  
  marker_plot_data <- marker_hr_table %>%
    arrange(HR) %>%
    mutate(Variable = factor(Variable, levels = unique(Variable)))
  
  marker_hr_plot <- ggplot(marker_plot_data, aes(x = Variable, y = HR)) +
    geom_point(color = sex_colors[["Male"]], size = 3) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper),
      width = 0.2,
      color = sex_colors[["Male"]]
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_y_log10() +
    labs(
      x = "Variables",
      y = "Hazard Ratio"
    ) +
    facet_wrap(~paste0(pred, " HRs with 90% CIs (", data_type, ")")) +
    theme_bw()
  
  save_pdf(
    marker_hr_plot,
    filename = paste0(pred, "_MaleDerived_", data_type, "_hr_plot.pdf"),
    width = 4.5,
    height = 3.5
  )
  
  # --------------------------------------------------------------------------
  # Adjusted survival plot using median-dichotomized marker
  # --------------------------------------------------------------------------
  
  new_data <- data.frame(
    cl = factor(c("High", "Low"), levels = c("Low", "High"))
  )
  
  rownames(new_data) <- paste0(
    c("High", "Low"),
    " N = ",
    table(data$cl)[c("High", "Low")]
  )
  
  adjustment_data <- make_adjustment_data(
    adjustment_covariates,
    nrow(new_data)
  )
  
  new_data <- data.frame(new_data, adjustment_data)
  
  marker_survfit <- survfit(
    fit_binary,
    newdata = new_data
  )
  
  continuous_p <- round(summary(fit_continuous)$coef[pred, 5], 3)
  binary_p <- round(summary(fit_binary)$coef["clHigh", 5], 3)
  
  splot <- ggsurvplot(
    marker_survfit,
    data = new_data,
    title = paste0(pred, " - ", data_type, " (Male-derived marker)"),
    conf.int = FALSE,
    pval = paste0(
      pred,
      " continuous pval: ",
      continuous_p
    ),
    palette = c(sex_colors[["Male"]], "gray40"),
    risk.table = TRUE,
    risk.table.fontsize = 5,
    xlab = "Time",
    ylab = "Survival Probability",
    legend.title = pred,
    legend.labs = c("High", "Low"),
    ggtheme = theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = text_size),
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size)
      ),
    tables.theme = theme_bw()
  )
  
  aligned <- cowplot::align_plots(
    splot$plot + theme(plot.margin = margin(r = 10)),
    splot$table + theme(plot.margin = margin(r = 10)),
    align = "v"
  )
  
  pdf(
    file.path(
      output_dir,
      paste0(pred, "_MaleDerived_", data_type, "_splot.pdf")
    ),
    width = 6,
    height = 7
  )
  
  print(
    cowplot::plot_grid(
      aligned[[1]],
      aligned[[2]],
      ncol = 1,
      rel_heights = c(4, 1)
    )
  )
  
  dev.off()
  
  invisible(
    list(
      marker = pred,
      fit_continuous = fit_continuous,
      fit_binary = fit_binary,
      anova_p = anova_p,
      hr_table = marker_hr_table
    )
  )
}

# ------------------------------------------------------------------------------
# Run MORC2 and PHKG2 validation analyses
# ------------------------------------------------------------------------------

single_marker_validation_results <- lapply(
  marker_genes,
  plot_single_marker_validation,
  data = validation_data,
  adjustment_covariates = adjustment_covariates
)

names(single_marker_validation_results) <- marker_genes

save(
  validation_data,
  single_marker_validation_results,
  file = file.path(
    output_dir,
    "single_marker_validation_results_MORC2_PHKG2.RData"
  )
)

message("Done. MORC2/PHKG2 validation outputs written to: ", output_dir)
