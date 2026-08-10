#!/usr/bin/env Rscript
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

# Reproduce the Figure 2C protein normal/tumor trajectory panel using temporalCPSA.
#
# Inputs:
# - Tumor protein matrix: ../data/cDisc_proteome_imputed_data_09152023.tsv
# - Clinical metadata: ../data/STable1.xlsx sheet ClinicalTable
# - Normal reference data: temporalCPSA package extdata normal_reference.rds
#
# Output:
# - Figure2C_protein_curveplot.pdf
# - Figure2C_protein_curveplot.png

if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA before running this script.", call. = FALSE)
}

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

data_dir <- "../data"
genes <- c("CNTN1", "MAPT", "L1CAM")
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
script_dir <- resolve_script_dir("Figure2C_protein_curveplot.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, "Figure2C_protein_curveplot.pdf")
output_png <- file.path(output_dir, "Figure2C_protein_curveplot.png")

# Figure 2C used the adaptive-span protein_tadj50 run in the original code.
# These are the feature/sex/tissue-specific spans stored in
# protein_tadj50_list.RData for CNTN1, MAPT, and L1CAM. Keeping these spans
# explicit is necessary for numerical equivalence with tn_df$value/low/hi.
figure2c_spans <- data.frame(
  feature = rep(genes, each = 4),
  sex = rep(c("Male", "Male", "Female", "Female"), times = length(genes)),
  tissue = rep(c("Normal", "Tumor", "Normal", "Tumor"), times = length(genes)),
  span = c(
    1.0, 1.3, 1.0, 0.5,
    1.1, 1.7, 1.3, 1.3,
    1.0, 1.3, 1.4, 0.5
  ),
  stringsAsFactors = FALSE
)

message("Loading clinical data...")
clinical_raw <- temporalCPSA::ageTMP_load_clinical(data_dir)
clinical <- data.frame(
  id = temporalCPSA::ageTMP_normalize_sample_ids(clinical_raw$id),
  age = clinical_raw$cDisc_age,
  sex = clinical_raw$cDisc_Gender,
  age_class = clinical_raw$cDisc_age_class_name_derived,
  stringsAsFactors = FALSE
)

message("Loading and collapsing public tumor protein data...")
protein_raw <- temporalCPSA::ageTMP_load_molecular(data_dir, modality = "protein")
protein <- temporalCPSA::ageTMP_split_annotation_matrix(
  protein_raw,
  annotation_cols = 1:4,
  row_id = "ApprovedGeneSymbol"
)
tumor_protein <- temporalCPSA::ageTMP_collapse_matrix_by_feature(
  protein$matrix,
  protein$annotation$ApprovedGeneSymbol
)

# Match proteo_tadj50.R: fit trajectories using tumor samples up to
# age.cut2 = 80. The public data include one extra male sample, C3N-02255
# at age 88, which was not part of the original protein_tadj50_list.RData fit.
fit_ids <- clinical$id[!is.na(clinical$age) & clinical$age <= 80]
tumor_protein <- tumor_protein[, intersect(colnames(tumor_protein), fit_ids), drop = FALSE]
clinical <- clinical[clinical$id %in% colnames(tumor_protein), , drop = FALSE]

message("Loading package normal reference data...")
normal_reference <- temporalCPSA::ageTMP_load_normal_reference()
normal_protein <- normal_reference$protein$matrix
normal_metadata <- normal_reference$protein$sample_metadata

# Match get_tn_df.R: tn.df curves are evaluated at all tumor sample ages
# <= 50. The male and female models are sex-stratified, but both are predicted
# on this common age support.
prediction_ages <- sort(clinical$age[
  clinical$id %in% colnames(tumor_protein) &
    !is.na(clinical$age) &
    clinical$age <= 50
])

message("Fitting normal/tumor trajectory comparison...")
trajectory <- temporalCPSA::ageTMP_compare_normal_tumor_trajectory(
  tumor_mat = tumor_protein,
  tumor_metadata = clinical,
  normal_mat = normal_protein,
  normal_metadata = normal_metadata,
  features = genes,
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  normal_sample_col = "ID",
  normal_age_col = "Age",
  normal_sex_col = "Gender",
  normal_covariates = c("pH", "PMI", "Ethnicity"),
  center_age_range = c(0, 50),
  span = figure2c_spans,
  prediction_ages = prediction_ages
)

trajectory <- trajectory[!is.na(trajectory$fit), , drop = FALSE]
trajectory$feature <- factor(trajectory$feature, levels = genes)
trajectory$sex <- factor(trajectory$sex, levels = c("Male", "Female"))
trajectory$tissue <- factor(trajectory$tissue, levels = c("Tumor", "Normal"))
trajectory$trajectory_group <- factor(
  paste(trajectory$sex, trajectory$tissue),
  levels = c("Male Tumor", "Male Normal", "Female Tumor", "Female Normal")
)

plot_lines <- trajectory[order(trajectory$feature, trajectory$sex, trajectory$tissue, trajectory$age), ]

# Match the tumor/normal sex-specific palette used in code_dump/tn_boxplot.R.
trajectory_colors <- c(
  "Male Tumor" = "#0707CF",
  "Male Normal" = "#A3A3FF",
  "Female Tumor" = "#CC0303",
  "Female Normal" = "#FFA7A7"
)

figure2c <- ggplot2::ggplot(
  plot_lines,
  ggplot2::aes(
    x = age,
    y = fit,
    group = interaction(trajectory_group, feature),
    color = trajectory_group,
    fill = trajectory_group
  )
) +
  ggplot2::geom_vline(
    xintercept = c(15, 26, 40),
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey50"
  ) +
  ggplot2::geom_line(linewidth = 3) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.5,
    color = NA
  ) +
  ggplot2::facet_grid(. ~ feature + sex, scales = "free", drop = TRUE) +
  ggplot2::coord_cartesian(xlim = c(0, 50), ylim = c(-2, 2)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  ggplot2::scale_color_manual(values = trajectory_colors, drop = FALSE) +
  ggplot2::scale_fill_manual(values = trajectory_colors, drop = FALSE) +
  ggplot2::labs(
    x = "age",
    y = "Tumor AD-MP",
    color = "Tissue",
    fill = "Tissue"
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 15),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    legend.title = ggplot2::element_text(size = 7),
    legend.text = ggplot2::element_text(size = 7),
    axis.text = ggplot2::element_text(size = 7),
    axis.title = ggplot2::element_text(size = 8),
    panel.spacing.x = grid::unit(1.8, "mm")
  )

ggplot2::ggsave(output_pdf, figure2c, width = 10, height = 3, units = "in")
ggplot2::ggsave(output_png, figure2c, width = 10, height = 3, units = "in", dpi = 220)

message("Saved ", output_pdf)
message("Saved ", output_png)
