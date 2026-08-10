#!/usr/bin/env Rscript
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# Purpose: Recreate the Figure 2F OXPHOS normal/tumor trajectory bundle from
#          public protein data using temporalCPSA.
# Input: data/cDisc_proteome_imputed_data_09152023.tsv,
#        data/STable1.xlsx sheet ClinicalTable,
#        temporalCPSA package normal reference data
# Output: Figure2F_oxphos_from_temporalCPSA.pdf and Figure2F_oxphos_from_temporalCPSA.png

if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA before running this script.", call. = FALSE)
}

required_packages <- c("ggplot2", "geomtextpath")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

data_dir <- "../data"
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
script_dir <- resolve_script_dir("Figure2F_oxphos_curveplot.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, "Figure2F_oxphos_from_temporalCPSA.pdf")
output_png <- file.path(output_dir, "Figure2F_oxphos_from_temporalCPSA.png")
apply_sd_filter <- FALSE
sd_min <- 0.15

# Figure 2F visualizes age-dependent tumor and normal/reference trajectories
# for proteins in the MITO3_OXPHOS module. This focused pathway definition is
# used instead of the broader HALLMARK_OXIDATIVE_PHOSPHORYLATION set so that
# the panel emphasizes the mitochondrial oxidative phosphorylation trajectory
# program. Trajectories are estimated below from the public tumor protein matrix
# and the package normal reference data.
oxphos_genes <- c(
  "ACAD9", "AIFM1", "ATP5F1A", "ATP5F1B", "ATP5F1C",
  "ATP5F1D", "ATP5F1E", "ATP5IF1", "ATP5MC1", "ATP5MC2",
  "ATP5MC3", "ATP5MK", "ATP5ME", "ATP5MF", "ATP5MG",
  "ATP5MJ", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
  "ATPAF1", "ATPAF2", "ATPSCKMT", "BCS1L", "CEP89",
  "CMC1", "CMC2", "COA1", "COA3", "COA4", "COA5", "COA6",
  "COA7", "COA8", "COX10", "COX11", "COX14", "COX15",
  "COX16", "COX17", "COX18", "COX19", "COX20", "COX4I1",
  "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2",
  "COX6B1", "COX6B2", "COX6C", "COX7A1", "COX7A2",
  "COX7A2L", "COX7B", "COX7B2", "COX7C", "COX8A",
  "COX8C", "CYC1", "CYCS", "DMAC1", "DMAC2", "DMAC2L",
  "ECSIT", "FMC1", "FOXRED1", "HCCS", "HIGD1A", "HIGD2A",
  "LYRM2", "LYRM7", "MT-ATP6", "MT-ATP8", "MT-CO1",
  "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2",
  "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
  "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13",
  "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6",
  "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFAF1",
  "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6",
  "NDUFAF7", "NDUFAF8", "NDUFB1", "NDUFB10", "NDUFB11",
  "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6",
  "NDUFB7", "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2",
  "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5",
  "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2",
  "NDUFV3", "NUBPL", "PET100", "PET117", "PNKD", "RAB5IF",
  "SCO1", "SCO2", "SDHA", "SDHAF1", "SDHAF2", "SDHAF3",
  "SDHAF4", "SDHB", "SDHC", "SDHD", "SMIM20", "SURF1",
  "TACO1", "TIMM21", "TIMMDC1", "TMEM126A", "TMEM126B",
  "TMEM177", "TMEM186", "TMEM70", "TTC19", "UQCC1",
  "UQCC2", "UQCC3", "UQCR10", "UQCR11", "UQCRB", "UQCRC1",
  "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ"
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

# Estimate trajectories on tumor samples within the analysis age range, then
# display the childhood-to-young-adult window emphasized in the OXPHOS panel.
fit_ids <- clinical$id[!is.na(clinical$age) & clinical$age <= 80]
tumor_protein <- tumor_protein[, intersect(colnames(tumor_protein), fit_ids), drop = FALSE]
clinical <- clinical[clinical$id %in% colnames(tumor_protein), , drop = FALSE]

message("Loading package normal reference data...")
normal_reference <- temporalCPSA::ageTMP_load_normal_reference()
normal_protein <- normal_reference$protein$matrix
normal_metadata <- normal_reference$protein$sample_metadata

genes <- intersect(oxphos_genes, intersect(rownames(tumor_protein), rownames(normal_protein)))
if (!"COX6C" %in% genes) {
  stop("COX6C is not present in both tumor and normal protein matrices.", call. = FALSE)
}
if (length(genes) < 6) {
  stop("Too few OXPHOS genes overlap the tumor and normal protein matrices.", call. = FALSE)
}
message("Fitting trajectories for ", length(genes), " OXPHOS proteins...")

# Predict trajectories at observed tumor-sample ages within the PED/ADO/YA
# support, including ages beyond the final displayed range. The ggplot loess
# layer then smooths the displayed trajectory bundle consistently across the
# visible age window.
prediction_ages <- sort(unique(clinical$age[
  !is.na(clinical$age) &
    !is.na(clinical$age_class) &
    clinical$age <= 49
]))

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
  adaptive_span = TRUE,
  tumor_min_span = 0.5,
  tumor_max_span = 3,
  normal_min_span = 1,
  normal_max_span = 3,
  span_step = 0.1,
  prediction_ages = prediction_ages
)

trajectory <- trajectory[!is.na(trajectory$fit), , drop = FALSE]
trajectory$sex <- factor(trajectory$sex, levels = c("Male", "Female"))
trajectory$tissue <- factor(trajectory$tissue, levels = c("Tumor", "Normal"))
trajectory$label <- factor(
  paste(trajectory$sex, trajectory$tissue),
  levels = c("Male Tumor", "Male Normal", "Female Tumor", "Female Normal")
)

trajectory_colors <- c(
  "Male Tumor" = "#0707CF",
  "Male Normal" = "#A3A3FF",
  "Female Tumor" = "#CC0303",
  "Female Normal" = "#FFA7A7"
)
trajectory_linetypes <- c(
  "Male Tumor" = "solid",
  "Male Normal" = "dashed",
  "Female Tumor" = "solid",
  "Female Normal" = "dashed"
)

sd_rank <- temporalCPSA::ageTMP_rank_trajectory_sd(
  trajectory,
  group_cols = "sex",
  tissue = "Tumor"
)
if (apply_sd_filter) {
  keep_genes <- temporalCPSA::ageTMP_filter_trajectory_sd(
    sd_rank,
    sd_min = sd_min,
    group_cols = "sex",
    keep = "all_groups"
  )
  trajectory <- trajectory[trajectory$feature %in% keep_genes, , drop = FALSE]
  message("Applied tumor trajectory SD filter sd >= ", sd_min, " in both sexes.")
  message("Retained ", length(unique(trajectory$feature)), " OXPHOS proteins after SD filtering.")
} else {
  message(
    "Not applying trajectory SD filtering; plotting all overlapping OXPHOS ",
    "proteins used for this trajectory bundle."
  )
}

label_data <- trajectory[
  trajectory$feature == "COX6C" &
    (
      (trajectory$tissue == "Tumor" & trajectory$age >= 30 & trajectory$age <= 37) |
        (trajectory$tissue == "Normal" & trajectory$age >= 28 & trajectory$age <= 36)
    ),
  ,
  drop = FALSE
]

figure2f <- ggplot2::ggplot(
  trajectory,
  ggplot2::aes(
    x = age,
    y = fit,
    group = interaction(feature, label),
    color = label,
    linetype = label
  )
) +
  ggplot2::geom_vline(
    xintercept = c(15, 26, 40),
    linetype = "dashed",
    linewidth = 0.25,
    color = "grey25"
  ) +
  ggplot2::geom_smooth(linewidth = 0.22, alpha = 0.95, se = FALSE) +
  geomtextpath::geom_labelsmooth(
    data = label_data,
    ggplot2::aes(label = feature),
    method = "loess",
    span = 1,
    se = FALSE,
    straight = FALSE,
    size = 3,
    linewidth = 0.35,
    fill = "white",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~sex, nrow = 1) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_ages), 40),
    ylim = c(-1.1, 1.1)
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(10, 20, 30, 40),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  ggplot2::scale_color_manual(values = trajectory_colors, drop = FALSE, name = "Tissue") +
  ggplot2::scale_linetype_manual(values = trajectory_linetypes, drop = FALSE, name = "Tissue") +
  ggplot2::labs(
    title = "OXPHOS",
    x = "age",
    y = "T-TMP or N-TMP value"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey85", color = "grey35", linewidth = 0.4),
    strip.text = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(size = 12, hjust = 0),
    axis.title = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(size = 9),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.width = grid::unit(1.1, "cm"),
    panel.spacing.x = grid::unit(2, "mm"),
    plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5)
  )

ggplot2::ggsave(output_pdf, figure2f, width = 8.6, height = 4.2, units = "in", useDingbats = FALSE)
ggplot2::ggsave(output_png, figure2f, width = 8.6, height = 4.2, units = "in", dpi = 300)

message("Saved ", output_pdf)
message("Saved ", output_png)
