#!/usr/bin/env Rscript
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# Purpose: Recreate the Figure 2G REACTOME_NEURONAL_SYSTEM normal/tumor trajectory bundle from
#          public protein data using temporalCPSA.
# Input: data/cDisc_proteome_imputed_data_09152023.tsv,
#        data/STable1.xlsx sheet ClinicalTable,
#        temporalCPSA package normal reference data
# Output: Figure2G_neuronal_system_from_temporalCPSA.pdf and Figure2G_neuronal_system_from_temporalCPSA.png

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
script_dir <- resolve_script_dir("Figure2G_neuronal_system_curveplot.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, "Figure2G_neuronal_system_from_temporalCPSA.pdf")
output_png <- file.path(output_dir, "Figure2G_neuronal_system_from_temporalCPSA.png")
apply_sd_filter <- FALSE
sd_min <- 0.15

# Figure 2G visualizes age-dependent tumor and normal/reference trajectories
# for proteins in the MSigDB v7.2 REACTOME_NEURONAL_SYSTEM gene set. This
# pathway definition captures neuronal signaling and synaptic programs that
# vary across developmental and tumor age trajectories. Trajectories are
# estimated below from the public tumor protein matrix and the package normal
# reference data.
neuronal_system_genes <- c(
  "ABAT", "ABCC8", "ABCC9", "ACHE", "ACTN2", "ADCY1", "ADCY2", "ADCY3", "ADCY4", "ADCY5",
  "ADCY6", "ADCY7", "ADCY8", "ADCY9", "AKAP5", "ALDH2", "ALDH5A1", "AP2A1", "AP2A2", "AP2B1",
  "AP2M1", "AP2S1", "APBA1", "APBA2", "APBA3", "ARHGEF7", "ARHGEF9", "ARL6IP5", "BCHE", "BEGAIN",
  "CACNA1A", "CACNA1B", "CACNA1E", "CACNA2D2", "CACNA2D3", "CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG2",
  "CACNG3", "CACNG4", "CACNG8", "CALM1", "CAMK1", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4",
  "CAMKK1", "CAMKK2", "CASK", "CHAT", "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6",
  "CHRNA7", "CHRNA9", "CHRNB2", "CHRNB3", "CHRNB4", "CHRND", "CHRNE", "CHRNG", "COMT", "CPLX1",
  "CREB1", "DBNL", "DLG1", "DLG2", "DLG3", "DLG4", "DLGAP1", "DLGAP2", "DLGAP3", "DLGAP4",
  "DNAJC5", "EPB41", "EPB41L1", "EPB41L2", "EPB41L3", "EPB41L5", "ERBB4", "FLOT1", "FLOT2", "GABBR1",
  "GABBR2", "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3",
  "GABRG2", "GABRG3", "GABRQ", "GABRR1", "GABRR2", "GABRR3", "GAD1", "GAD2", "GIT1", "GJA10",
  "GJC1", "GJD2", "GLRA1", "GLRA2", "GLRA3", "GLRB", "GLS", "GLS2", "GLUL", "GNAI1",
  "GNAI2", "GNAI3", "GNAL", "GNAT3", "GNB1", "GNB2", "GNB3", "GNB4", "GNB5", "GNG10",
  "GNG11", "GNG12", "GNG13", "GNG2", "GNG3", "GNG4", "GNG5", "GNG7", "GNG8", "GNGT1",
  "GNGT2", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
  "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "GRIP1", "GRIP2", "GRM1",
  "GRM5", "HCN1", "HCN2", "HCN3", "HCN4", "HOMER1", "HOMER2", "HOMER3", "HRAS", "HSPA8",
  "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E", "IL1RAP", "IL1RAPL1", "IL1RAPL2", "KCNA1", "KCNA10",
  "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6", "KCNA7", "KCNAB1", "KCNAB2", "KCNAB3", "KCNB1",
  "KCNB2", "KCNC1", "KCNC2", "KCNC3", "KCNC4", "KCND1", "KCND2", "KCND3", "KCNF1", "KCNG1",
  "KCNG2", "KCNG3", "KCNG4", "KCNH1", "KCNH2", "KCNH3", "KCNH4", "KCNH5", "KCNH6", "KCNH7",
  "KCNH8", "KCNJ1", "KCNJ10", "KCNJ11", "KCNJ12", "KCNJ14", "KCNJ15", "KCNJ16", "KCNJ2", "KCNJ3",
  "KCNJ4", "KCNJ5", "KCNJ6", "KCNJ8", "KCNJ9", "KCNK1", "KCNK10", "KCNK13", "KCNK16", "KCNK17",
  "KCNK18", "KCNK2", "KCNK3", "KCNK4", "KCNK6", "KCNK7", "KCNK9", "KCNMA1", "KCNMB1", "KCNMB2",
  "KCNMB3", "KCNMB4", "KCNN1", "KCNN2", "KCNN3", "KCNN4", "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4",
  "KCNQ5", "KCNS1", "KCNS2", "KCNS3", "KCNV1", "KCNV2", "KIF17", "KPNA2", "KRAS", "LIN7A",
  "LIN7B", "LIN7C", "LRFN1", "LRFN2", "LRFN3", "LRFN4", "LRRC4B", "LRRC7", "LRRTM1", "LRRTM2",
  "LRRTM3", "LRRTM4", "LRTOMT", "MAOA", "MAPK1", "MAPK3", "MAPT", "MDM2", "MYO6", "NAAA",
  "NBEA", "NCALD", "NEFL", "NLGN1", "NLGN2", "NLGN3", "NLGN4X", "NLGN4Y", "NPTN", "NRAS",
  "NRG1", "NRGN", "NRXN1", "NRXN2", "NRXN3", "NSF", "NTRK3", "PANX1", "PANX2", "PDLIM5",
  "PDPK1", "PICK1", "PLCB1", "PLCB2", "PLCB3", "PPFIA1", "PPFIA2", "PPFIA3", "PPFIA4", "PPFIBP1",
  "PPFIBP2", "PPM1E", "PPM1F", "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKACA", "PRKACB", "PRKACG",
  "PRKAG1", "PRKAG2", "PRKAG3", "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKCA", "PRKCB", "PRKCG",
  "PRKX", "PTPRD", "PTPRF", "PTPRS", "RAB3A", "RAC1", "RASGRF1", "RASGRF2", "RIMS1", "RPS6KA1",
  "RPS6KA2", "RPS6KA3", "RPS6KA6", "RTN3", "SHANK1", "SHANK2", "SHARPIN", "SIPA1L1", "SLC17A7", "SLC18A2",
  "SLC18A3", "SLC1A1", "SLC1A2", "SLC1A3", "SLC1A6", "SLC1A7", "SLC22A1", "SLC22A2", "SLC32A1", "SLC38A1",
  "SLC38A2", "SLC5A7", "SLC6A1", "SLC6A11", "SLC6A12", "SLC6A13", "SLC6A3", "SLC6A4", "SLITRK1", "SLITRK2",
  "SLITRK3", "SLITRK4", "SLITRK5", "SLITRK6", "SNAP25", "SRC", "STX1A", "STXBP1", "SYN1", "SYN2",
  "SYN3", "SYT1", "SYT10", "SYT12", "SYT2", "SYT7", "SYT9", "TSPAN7", "TSPOAP1", "TUBA1A",
  "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D", "TUBA3E", "TUBA4A", "TUBA4B", "TUBA8", "TUBAL3", "TUBB1",
  "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B", "TUBB6", "TUBB8", "TUBB8B", "UNC13B", "VAMP2"
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
# display the childhood-to-young-adult window emphasized in the neuronal-system
# panel.
fit_ids <- clinical$id[!is.na(clinical$age) & clinical$age <= 80]
tumor_protein <- tumor_protein[, intersect(colnames(tumor_protein), fit_ids), drop = FALSE]
clinical <- clinical[clinical$id %in% colnames(tumor_protein), , drop = FALSE]

message("Loading package normal reference data...")
normal_reference <- temporalCPSA::ageTMP_load_normal_reference()
normal_protein <- normal_reference$protein$matrix
normal_metadata <- normal_reference$protein$sample_metadata

genes <- intersect(neuronal_system_genes, intersect(rownames(tumor_protein), rownames(normal_protein)))
if (!"PRKAR2B" %in% genes) {
  stop("PRKAR2B is not present in both tumor and normal protein matrices.", call. = FALSE)
}
if (length(genes) < 6) {
  stop("Too few neuronal-system genes overlap the tumor and normal protein matrices.", call. = FALSE)
}
message("Fitting trajectories for ", length(genes), " neuronal-system proteins...")

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
  message("Retained ", length(unique(trajectory$feature)), " neuronal-system proteins after SD filtering.")
} else {
  message(
    "Not applying trajectory SD filtering; plotting all overlapping neuronal-system ",
    "proteins used for this trajectory bundle."
  )
}

highlight_genes <- c("PRKAR2B", "PRKACB")
label_data <- trajectory[
  (
    trajectory$feature == "PRKAR2B" & trajectory$age > 28 & trajectory$age < 38
  ) |
    (
      trajectory$feature == "PRKACB" & trajectory$age > 20 & trajectory$age < 25
    ),
  ,
  drop = FALSE
]

figure2g <- ggplot2::ggplot(
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
  ggplot2::geom_smooth(
    data = trajectory[trajectory$feature %in% highlight_genes, , drop = FALSE],
    linewidth = 0.7,
    alpha = 1,
    se = FALSE,
    show.legend = FALSE
  ) +
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
    title = "REACTOME_NEURONAL_SYSTEM",
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

ggplot2::ggsave(output_pdf, figure2g, width = 8.6, height = 4.2, units = "in", useDingbats = FALSE)
ggplot2::ggsave(output_png, figure2g, width = 8.6, height = 4.2, units = "in", dpi = 300)

message("Saved ", output_pdf)
message("Saved ", output_png)
