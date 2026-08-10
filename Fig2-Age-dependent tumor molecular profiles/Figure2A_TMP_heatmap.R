#!/usr/bin/env Rscript
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# Purpose: Recreate Figure 2A tumor/normal AD-TMP heatmap from public data using temporalCPSA.
#
# Design note:
# Figure 2A is a coordinated tumor/reference trajectory heatmap. Its row order
# encodes the relationship among sex-specific tumor trajectory groups and
# normal/reference trajectory groups. The T-TMP and N-TMP matrices remain
# gene-aligned row by row, while genes are ordered by a four-part alignment key:
# male tumor group, male normal group, female tumor group, and female normal
# group. This makes the sex-specific tumor/reference trajectory structure
# visible without relying on serialized intermediate objects.
# Input: data/STable1.xlsx sheet ClinicalTable, data/cDisc_proteome_imputed_data_09152023.tsv,
#        data/STable2.xlsx sheet TTMP_Protein_Group.
# Output: Figure2A_from_temporalCPSA.pdf and Figure2A_from_temporalCPSA.png.

required <- c("ComplexHeatmap", "circlize", "grid", "readxl", "RColorBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA before running this script.", call. = FALSE)
}

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
script_dir <- resolve_script_dir("Figure2A_TMP_heatmap.R")
`%||%` <- function(x, y) if (is.null(x)) y else x

repo_data <- Sys.getenv("FIGURE_DATA_DIR", file.path(script_dir, "..", "data"))
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
clinical_file <- file.path(repo_data, "STable1.xlsx")
protein_file <- file.path(repo_data, "cDisc_proteome_imputed_data_09152023.tsv")
stable2_file <- file.path(repo_data, "STable2.xlsx")
for (file in c(clinical_file, protein_file, stable2_file)) {
  if (!file.exists(file)) {
    stop("Required input file not found: ", file, call. = FALSE)
  }
}

stable2_sheets <- readxl::excel_sheets(stable2_file)
stable1_sheets <- readxl::excel_sheets(clinical_file)
if (!"ClinicalTable" %in% stable1_sheets) {
  stop("Required STable1 sheet not found: ClinicalTable", call. = FALSE)
}
for (sheet in "TTMP_Protein_Group") {
  if (!sheet %in% stable2_sheets) {
    stop("Required STable2 sheet not found: ", sheet, call. = FALSE)
  }
}
cl2_sankey_sheet_candidates <- c(
  "cl2_sankey_data",
  "cl2.sankey.data",
  "cl2 sankey data",
  "cl2_snakey_data"
)
cl2_sankey_sheet <- stable2_sheets[
  match(tolower(cl2_sankey_sheet_candidates), tolower(stable2_sheets), nomatch = 0)
]
cl2_sankey_sheet <- cl2_sankey_sheet[cl2_sankey_sheet != ""][1] %||% NA_character_

standardize_cl2_sankey_columns <- function(x) {
  original_names <- names(x)
  normalized_names <- tolower(gsub("[^a-z0-9]+", "_", original_names))
  normalized_names <- gsub("^_|_$", "", normalized_names)

  rename_one <- function(target, aliases) {
    hit <- which(normalized_names %in% aliases)
    if (length(hit) > 0 && !target %in% names(x)) {
      names(x)[hit[1]] <<- target
    }
  }

  rename_one("gene", c("gene", "genes", "protein", "protein_gene", "symbol"))
  rename_one("mf.cl", c("mf_cl", "male_female_cl", "male_female_t_tmp_group"))
  rename_one("male.cl", c("male_cl", "male_t_tmp_group", "male_group"))
  rename_one("female.cl", c("female_cl", "female_t_tmp_group", "female_group"))
  rename_one("n.male.cl", c("n_male_cl", "male_n_tmp_group", "male_normal_group", "normal_male_cl"))
  rename_one("n.female.cl", c("n_female_cl", "female_n_tmp_group", "female_normal_group", "normal_female_cl"))

  x
}

message("Loading public clinical, protein, and STable2 data...")
clinical <- temporalCPSA::ageTMP_load_clinical(repo_data)
protein_raw <- temporalCPSA::ageTMP_load_molecular(repo_data, "protein")
groups <- as.data.frame(temporalCPSA::ageTMP_load_supplement(repo_data, "STable2", "TTMP_Protein_Group"))
groups <- standardize_cl2_sankey_columns(groups)
if ("gene" %in% names(groups) && !"Gene" %in% names(groups)) {
  names(groups)[names(groups) == "gene"] <- "Gene"
}
if ("male.cl" %in% names(groups) && !"Male T-TMP Group" %in% names(groups)) {
  groups[["Male T-TMP Group"]] <- groups[["male.cl"]]
}
if ("female.cl" %in% names(groups) && !"Female T-TMP Group" %in% names(groups)) {
  groups[["Female T-TMP Group"]] <- groups[["female.cl"]]
}
if ("mf.cl" %in% names(groups) && !"Male:Female T-TMP Group" %in% names(groups)) {
  groups[["Male:Female T-TMP Group"]] <- groups[["mf.cl"]]
}

required_clinical <- c("id", "cDisc_age", "cDisc_Gender", "cDisc_age_class_name_derived")
missing_clinical <- setdiff(required_clinical, names(clinical))
if (length(missing_clinical) > 0) {
  stop("Clinical data missing required column(s): ", paste(missing_clinical, collapse = ", "), call. = FALSE)
}
required_group <- c("Gene", "Male T-TMP Group", "Female T-TMP Group", "Male:Female T-TMP Group")
missing_group <- setdiff(required_group, names(groups))
if (length(missing_group) > 0) {
  stop("TTMP_Protein_Group missing required column(s): ", paste(missing_group, collapse = ", "), call. = FALSE)
}
required_protein <- c("ApprovedGeneSymbol", "Symbol.V5")
missing_protein <- setdiff(required_protein, names(protein_raw))
if (length(missing_protein) > 0) {
  stop("Protein data missing required annotation column(s): ", paste(missing_protein, collapse = ", "), call. = FALSE)
}

message("Preparing gene-level protein matrix...")
annotation_cols <- match(c("ApprovedGeneSymbol", "Symbol.V5", "OldSymbol", "coding_gene"), names(protein_raw))
annotation_cols <- annotation_cols[!is.na(annotation_cols)]
protein_parts <- temporalCPSA::ageTMP_split_annotation_matrix(protein_raw, annotation_cols = annotation_cols)
protein_mat <- temporalCPSA::ageTMP_collapse_matrix_by_feature(
  protein_parts$matrix,
  protein_parts$annotation$ApprovedGeneSymbol
)

clinical$id <- temporalCPSA::ageTMP_normalize_sample_ids(clinical$id)
clinical$age <- as.numeric(clinical$cDisc_age)
clinical$sex <- clinical$cDisc_Gender
clinical$age_class <- clinical$cDisc_age_class_name_derived

sample_ids <- intersect(colnames(protein_mat), clinical$id)
protein_mat <- protein_mat[, sample_ids, drop = FALSE]
clinical <- clinical[match(sample_ids, clinical$id), , drop = FALSE]

# Protein AD-TMP trajectories are estimated on the analysis age range used for
# the Figure 2 trajectory panels. Restricting to age <= 80 keeps the heatmap on
# the same age support as the trajectory model used for this figure.
fit_keep <- !is.na(clinical$age) & clinical$age <= 80
clinical <- clinical[fit_keep, , drop = FALSE]
protein_mat <- protein_mat[, clinical$id, drop = FALSE]

age_levels <- c("PED", "ADO", "YA", "ADULT")
plot_clinical <- clinical[
  !is.na(clinical$age) &
    !is.na(clinical$age_class) &
    clinical$age_class %in% age_levels &
    clinical$age <= 47 &
    clinical$id %in% colnames(protein_mat),
  , drop = FALSE
]
plot_clinical$age_class <- factor(plot_clinical$age_class, levels = age_levels)
plot_clinical <- plot_clinical[order(plot_clinical$age), , drop = FALSE]

features <- intersect(groups$Gene, rownames(protein_mat))
groups <- groups[match(features, groups$Gene), , drop = FALSE]
if (length(features) == 0) {
  stop("No STable2 T-TMP genes overlap the public protein matrix.", call. = FALSE)
}
message("Restricting Figure 2A heatmap to proteins present in tumor and normal/reference matrices...")
normal_reference <- temporalCPSA::ageTMP_load_normal_reference()
normal_protein <- normal_reference$protein$matrix
normal_metadata <- normal_reference$protein$sample_metadata
groups_full <- groups
features_all_tumor <- intersect(groups_full$Gene, rownames(protein_mat))
features <- intersect(features_all_tumor, rownames(normal_protein))
groups <- groups_full[match(features, groups_full$Gene), , drop = FALSE]
if (length(features) == 0) {
  stop("No STable2 T-TMP genes overlap both the public protein matrix and package normal reference.", call. = FALSE)
}
message("Fitting tumor/normal AD-TMP trajectories for ", length(features), " proteins.")
message("Using adaptive feature/sex/tissue-specific loess spans for trajectory estimation.")

prediction_ages <- plot_clinical$age
prediction_ids <- plot_clinical$id
col_age_class <- plot_clinical$age_class

trajectory <- temporalCPSA::ageTMP_compare_normal_tumor_trajectory(
  tumor_mat = protein_mat,
  tumor_metadata = clinical,
  normal_mat = normal_protein,
  normal_metadata = normal_metadata,
  features = features,
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
  prediction_ages = prediction_ages,
  prediction_sample_ids = prediction_ids,
  ci_level = 0.95
)

finite_trajectory_fits <- sum(is.finite(trajectory$fit))
message("Finite trajectory fits: ", finite_trajectory_fits, " / ", length(trajectory$fit))
if (finite_trajectory_fits == 0) {
  stop("Trajectory fitting produced no finite AD-TMP values.", call. = FALSE)
}

make_fit_matrix <- function(traj, sex, tissue, feature_order, sample_order) {
  sex_df <- traj[traj$sex == sex & traj$tissue == tissue, c("feature", "sample_id", "fit"), drop = FALSE]
  out <- matrix(NA_real_, nrow = length(feature_order), ncol = length(sample_order),
                dimnames = list(feature_order, sample_order))
  split_df <- split(sex_df, sex_df$feature)
  for (feature in names(split_df)) {
    idx <- match(split_df[[feature]]$sample_id, sample_order)
    out[feature, idx] <- split_df[[feature]]$fit
  }
  out
}

cluster_normal_rows <- function(mat, k = 4) {
  scaled <- t(scale(t(mat)))
  scaled[!is.finite(scaled)] <- 0
  set.seed(1)
  km <- stats::kmeans(scaled, centers = k, nstart = 50, iter.max = 100)
  centers <- km$centers
  age_score <- rowMeans(centers[, seq_len(max(1, floor(ncol(centers) / 4))), drop = FALSE]) -
    rowMeans(centers[, seq.int(ceiling(ncol(centers) * 3 / 4), ncol(centers)), drop = FALSE])
  remap <- setNames(seq_len(k), names(sort(age_score, decreasing = TRUE)))
  out <- factor(remap[as.character(km$cluster)], levels = as.character(seq_len(k)))
  names(out) <- rownames(mat)
  out
}

message("Assigning male and female N-TMP trajectory groups...")
normal_cluster_feature_order <- groups$Gene
normal_consensus_k4 <- normal_reference$protein$clusters$consensus_k4 %||% NULL
if (!is.null(normal_consensus_k4)) {
  message("Using packaged N-TMP consensus clusters with predefined cluster remapping.")
  normal_cluster_by_sex <- list(
    Male = setNames(
      factor(
        as.character(normal_consensus_k4[match(paste0("Male,", normal_cluster_feature_order), names(normal_consensus_k4))]),
        levels = as.character(1:4)
      ),
      normal_cluster_feature_order
    ),
    Female = setNames(
      factor(
        as.character(normal_consensus_k4[match(paste0("Female,", normal_cluster_feature_order), names(normal_consensus_k4))]),
        levels = as.character(1:4)
      ),
      normal_cluster_feature_order
    )
  )
} else {
  message("Packaged N-TMP consensus clusters not found; estimating joint k-means clusters as fallback.")
  # The fallback clusters are estimated in one shared male+female space. This is
  # the computational counterpart of the row-alignment logic: a color such as
  # normal group 2 should have the same trajectory meaning in the male and
  # female normal matrices, rather than being independently discovered and
  # independently colored in each sex.
  male_normal_for_cluster <- make_fit_matrix(
    trajectory, "Male", "Normal", normal_cluster_feature_order, prediction_ids
  )
  female_normal_for_cluster <- make_fit_matrix(
    trajectory, "Female", "Normal", normal_cluster_feature_order, prediction_ids
  )
  joint_normal_for_cluster <- rbind(male_normal_for_cluster, female_normal_for_cluster)
  rownames(joint_normal_for_cluster) <- c(
    paste0("Male,", normal_cluster_feature_order),
    paste0("Female,", normal_cluster_feature_order)
  )
  joint_normal_clusters <- cluster_normal_rows(joint_normal_for_cluster, k = 4)
  normal_cluster_by_sex <- list(
    Male = setNames(
      joint_normal_clusters[match(paste0("Male,", normal_cluster_feature_order), names(joint_normal_clusters))],
      normal_cluster_feature_order
    ),
    Female = setNames(
      joint_normal_clusters[match(paste0("Female,", normal_cluster_feature_order), names(joint_normal_clusters))],
      normal_cluster_feature_order
    )
  )
}
normal_alignment_summary <- table(
  male = normal_cluster_by_sex$Male[normal_cluster_feature_order],
  female = normal_cluster_by_sex$Female[normal_cluster_feature_order]
)
print(normal_alignment_summary)

normal_cluster_lookup <- factor(
  c(
    as.character(normal_cluster_by_sex$Male[normal_cluster_feature_order]),
    as.character(normal_cluster_by_sex$Female[normal_cluster_feature_order])
  ),
  levels = c("3", "4", "1", "2")
)
names(normal_cluster_lookup) <- c(
  paste0("Male,", normal_cluster_feature_order),
  paste0("Female,", normal_cluster_feature_order)
)

get_green <- function(n) grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-(1:2)])(n)
cluster_col <- c("1" = "#ef1c1c", "2" = "#3a86c8", "3" = "#2db84d", "4" = "#7a00cc")
normal_cluster_col <- c("1" = "#F6CFC7", "2" = "#E686E9", "3" = "#EAA1B8", "4" = "#6A80D8")
age_class_col <- setNames(get_green(4), age_levels)
tissue_col <- c(Tumor = "#FFA500", Normal = "#A020F0")
sex_col <- c(Male = "#0000CC", Female = "#CC0000")
heatmap_col <- circlize::colorRamp2(c(-1, 0, 1), c("#2B00FF", "white", "#FF2A1A"))

# Curated pathway membership columns used in the long Figure 2A layout. These
# are visualization annotations; the trajectory matrices and cluster rows above
# are regenerated from public data and temporalCPSA.
plot_pathways <- c(
  "GOMF_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "MITO3_OXPHOS",
  "GOMF_CALMODULIN_BINDING",
  "GOBP_NEURON_DEVELOPMENT",
  "GOBP_CELL_CELL_SIGNALING",
  "REACTOME_NEURONAL_SYSTEM"
)
pathway_gene_sets <- list(
  GOMF_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY = c(
    "ATP1A3", "ATP2B2", "ATP6V0A1", "ATP6V0D1", "ATP6V1A",
    "ATP6V1B2", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F",
    "ATP6V1H", "CACNA2D1", "COX4I1", "COX5A", "COX5B", "COX6B1",
    "CYC1", "GPM6A", "NNT", "SFXN1", "SLC12A5", "SLC17A7",
    "SLC1A2", "SLC1A3", "SLC25A12", "SLC25A22", "SLC25A3",
    "SLC3A2", "SLC4A4", "SLC8A2", "SNAP25", "UQCRC1",
    "UQCRFS1", "UQCRH", "VDAC1", "VDAC2"
  ),
  MITO3_OXPHOS = c(
    "ABAT", "ACAT1", "ACO2", "ACOT7", "AFG3L2", "ALDH2",
    "ALDH4A1", "ALDH5A1", "ALDH6A1", "ALDH9A1", "APOO", "C1QBP",
    "CISD1", "COX4I1", "COX5A", "COX5B", "COX6B1", "COX6C",
    "COX7A2", "COX7C", "CS", "CYC1", "CYCS", "DCXR", "DLAT",
    "DLD", "DLST", "DNM1L", "ECH1", "FH", "FIS1", "GLS",
    "GOT2", "GPD2", "HADH", "HINT1", "HSPA9", "IDH3A", "IDH3B",
    "IMMT", "LDHB", "MDH2", "ME3", "MGST3", "NDUFA10",
    "NDUFA13", "NDUFA2", "NDUFA4", "NDUFA5", "NDUFA6",
    "NDUFA9", "NDUFAB1", "NDUFB10", "NDUFB4", "NDUFB6",
    "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5",
    "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2",
    "NIPSNAP1", "NNT", "OGDH", "OPA1", "OXCT1", "PDE2A",
    "PDHA1", "PDHB", "PHB1", "PHB2", "PRDX6", "PRXL2A",
    "QDPR", "SDHA", "SFXN1", "SLC25A11", "SLC25A12",
    "SLC25A22", "SLC25A3", "SLC25A6", "SOD2", "SUCLA2",
    "SUCLG1", "TOMM70", "UQCRC1", "UQCRC2", "UQCRFS1",
    "UQCRH", "UQCRQ", "VDAC1", "VDAC2"
  ),
  GOMF_CALMODULIN_BINDING = c(
    "ADD1", "ATP2B2", "CAMK2A", "CAMKV", "GAP43", "MAP6",
    "MBP", "MYO5A", "NRGN", "PLCB1", "PPP3CA", "PPP3R1",
    "SLC8A2", "SYT1", "VAMP2"
  ),
  GOBP_NEURON_DEVELOPMENT = c(
    "AFG3L2", "AP2A1", "ATL1", "ATP1B2", "CAMK2A", "CNP",
    "CNTN1", "CNTNAP1", "CPNE6", "CTNND2", "CYFIP2", "DCLK1",
    "DLG4", "EPB41L3", "GAP43", "GPM6A", "GPM6B", "L1CAM",
    "MAG", "MAP6", "MAPT", "MT3", "NCAM1", "NCAM2", "NCDN",
    "NDRG4", "NEGR1", "NPTN", "OGDH", "OMG", "OPCML", "PACSIN1",
    "PAK1", "PLP1", "PLXNA1", "PPP3CA", "RAB3A", "RAP2A",
    "SLC12A5", "SLC1A3", "SNAP25", "SOD2", "STX1B", "STXBP1",
    "SULT4A1", "SYN1", "SYT1", "THY1", "TNR", "UQCRQ", "VAPA",
    "WASF1"
  ),
  GOBP_CELL_CELL_SIGNALING = c(
    "ABAT", "CACNA2D1", "CALB2", "CAMK2A", "CAMKV", "CNP",
    "CPLX2", "CTNND2", "DLG4", "DNAJC5", "DNM1", "GLS", "GNAQ",
    "GNAZ", "HADH", "MAPT", "MBP", "NAPB", "NCDN", "NPTN",
    "PHF24", "PIN1", "PLCB1", "PLP1", "PPP3CA", "PPP3R1",
    "PRKCB", "PRKCG", "PRNP", "PRRT2", "RAB3A", "RAP1B",
    "SEPTIN5", "SLC12A5", "SLC17A7", "SLC1A2", "SLC1A3",
    "SLC25A22", "SLC8A2", "SNAP25", "STX1B", "STXBP1", "SV2A",
    "SV2B", "SYN1", "SYN3", "SYNGR1", "SYP", "SYT1", "THY1",
    "TMOD2", "TNR", "VAMP2", "VCP", "VDAC1", "VSNL1"
  ),
  REACTOME_NEURONAL_SYSTEM = c(
    "ABAT", "ALDH2", "ALDH5A1", "AP2A1", "AP2M1", "AP2S1",
    "CAMK2A", "DLG4", "DNAJC5", "EPB41L3", "GLS", "GLUL",
    "GNG2", "GNG3", "GNG7", "HSPA8", "MAPT", "NPTN", "NRGN",
    "NSF", "PLCB1", "PRKCB", "PRKCG", "RAB3A", "RTN3",
    "SLC17A7", "SLC1A2", "SLC1A3", "SNAP25", "STXBP1", "SYN1",
    "SYN3", "SYT1", "VAMP2"
  )
)

message("Assembling aligned normal/tumor Figure 2A heatmap...")
male_groups <- factor(as.character(groups[["Male T-TMP Group"]]), levels = as.character(1:4))
female_groups <- factor(as.character(groups[["Female T-TMP Group"]]), levels = as.character(1:4))

groups_has_cl2_sankey_data <- all(c(
  "Gene", "mf.cl", "male.cl", "female.cl", "n.male.cl", "n.female.cl"
) %in% names(groups))

if (!is.na(cl2_sankey_sheet) && cl2_sankey_sheet %in% stable2_sheets) {
  message("Using supplied STable2 ", cl2_sankey_sheet, " sheet for Figure 2A row alignment.")
  cl2_sankey_data <- standardize_cl2_sankey_columns(
    as.data.frame(readxl::read_excel(stable2_file, sheet = cl2_sankey_sheet))
  )
} else if (groups_has_cl2_sankey_data) {
  message("Using STable2 TTMP_Protein_Group full alignment data for Figure 2A row alignment.")
  cl2_sankey_data <- groups
  names(cl2_sankey_data)[names(cl2_sankey_data) == "Gene"] <- "gene"
} else {
  cl2_sankey_data <- NULL
}

if (!is.null(cl2_sankey_data)) {
  required_cl2 <- c(
    "gene", "mf.cl", "male.cl", "female.cl",
    "n.male.cl", "n.female.cl"
  )
  missing_cl2 <- setdiff(required_cl2, names(cl2_sankey_data))
  if (length(missing_cl2) > 0) {
    stop(
      "Supplied Figure 2A alignment table missing required column(s): ",
      paste(missing_cl2, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(cl2_sankey_data$gene) > 0) {
    stop("Supplied Figure 2A alignment table contains duplicate gene rows.", call. = FALSE)
  }
  cl2_sankey_data <- cl2_sankey_data[cl2_sankey_data$gene %in% features, , drop = FALSE]
  if (nrow(cl2_sankey_data) == 0) {
    stop("Supplied Figure 2A alignment table has no genes overlapping the Figure 2A features.", call. = FALSE)
  }
  missing_cl2_genes <- setdiff(features, cl2_sankey_data$gene)
  if (length(missing_cl2_genes) > 0) {
    stop(
      "Supplied Figure 2A alignment table is missing Figure 2A feature(s): ",
      paste(head(missing_cl2_genes, 20), collapse = ", "),
      if (length(missing_cl2_genes) > 20) " ..." else "",
      call. = FALSE
    )
  }
  for (cluster_col_name in c("male.cl", "female.cl")) {
    cl2_sankey_data[[cluster_col_name]] <- factor(
      as.character(cl2_sankey_data[[cluster_col_name]]),
      levels = as.character(1:4)
    )
  }
  # Normal trajectory clusters are ordered as 3, 4, 1, 2 for the row alignment
  # key. This maps normal/reference trajectory patterns into the visual order
  # used for the coordinated tumor-normal heatmap while retaining the displayed
  # n.cl labels as 1, 2, 3, and 4.
  for (cluster_col_name in c("n.male.cl", "n.female.cl")) {
    cl2_sankey_data[[cluster_col_name]] <- factor(
      as.character(cl2_sankey_data[[cluster_col_name]]),
      levels = c("3", "4", "1", "2")
    )
  }
  cl2_sankey_data$mf.cl <- as.character(cl2_sankey_data$mf.cl)
} else {
  message("Full TTMP_Protein_Group alignment columns not found; reconstructing row alignment from public groups and packaged N-TMP clusters.")
  # Reconstruct the alignment table used to coordinate the male/female tumor
  # groups with matched male/female normal trajectory groups. The row order is
  # generated in two steps: first order genes by the four-part alignment key,
  # then regroup duplicated male/female heatmap rows by sex, tumor group, and
  # normal/reference group.
  cl_sankey_data <- data.frame(
    gene = groups$Gene,
    male.cl = male_groups,
    female.cl = female_groups,
    mf.cl = groups[["Male:Female T-TMP Group"]],
    stringsAsFactors = FALSE
  )
  cl_sankey_data$n.female.cl <- normal_cluster_lookup[match(
    paste0("Female,", cl_sankey_data$gene),
    names(normal_cluster_lookup)
  )]
  cl_sankey_data$n.male.cl <- normal_cluster_lookup[match(
    paste0("Male,", cl_sankey_data$gene),
    names(normal_cluster_lookup)
  )]
  cl_sankey_data$n.mf.cl <- as.character(apply(
    cl_sankey_data[, c("n.male.cl", "n.female.cl"), drop = FALSE],
    1,
    function(x) ifelse(any(is.na(x)), NA_character_, paste0(x[1], x[2]))
  ))
  cl_sankey_data$nt.f <- paste0(cl_sankey_data$n.female.cl, cl_sankey_data$female.cl)
  cl_sankey_data$nt.m <- paste0(cl_sankey_data$n.male.cl, cl_sankey_data$male.cl)
  cl2_sankey_data <- cl_sankey_data[!is.na(cl_sankey_data$n.male.cl), , drop = FALSE]
  cl2_sankey_data$join <- factor(paste0(
    cl2_sankey_data$male.cl,
    cl2_sankey_data$n.male.cl,
    cl2_sankey_data$female.cl,
    cl2_sankey_data$n.female.cl
  ))
  cl2_sankey_data <- cl2_sankey_data[order(cl2_sankey_data$join), , drop = FALSE]
}

fit_mats <- list(
  Male = list(
    Normal = make_fit_matrix(trajectory, "Male", "Normal", groups$Gene, prediction_ids),
    Tumor = make_fit_matrix(trajectory, "Male", "Tumor", groups$Gene, prediction_ids)
  ),
  Female = list(
    Normal = make_fit_matrix(trajectory, "Female", "Normal", groups$Gene, prediction_ids),
    Tumor = make_fit_matrix(trajectory, "Female", "Tumor", groups$Gene, prediction_ids)
  )
)

pathway_col <- rep(list(c("TRUE" = "black", "FALSE" = "white")), length(plot_pathways))
names(pathway_col) <- plot_pathways

mark_genes <- c(
  "EPB41L3", "L1CAM", "NDUFB4", "CACNA2D1", "ATL1", "CYC1",
  "OPCML", "SLC8A2", "GNAZ", "MAPT", "NCAM2", "NDUFA6",
  "VDAC1", "CNTN1", "DCLK1", "CALB2"
)

row_df <- data.frame(
  label = c(
    paste0("Male,", cl2_sankey_data$gene),
    paste0("Female,", cl2_sankey_data$gene)
  ),
  stringsAsFactors = FALSE
)
row_df$gene <- gsub("(Male|Female),", "", row_df$label)
row_df$sex <- factor(gsub(",.*", "", row_df$label), levels = c("Male", "Female"))
row_df$join <- cl2_sankey_data[match(row_df$gene, cl2_sankey_data$gene), ]$mf.cl
row_df$male.cl <- cl2_sankey_data[match(row_df$gene, cl2_sankey_data$gene), ]$male.cl
row_df$female.cl <- cl2_sankey_data[match(row_df$gene, cl2_sankey_data$gene), ]$female.cl
row_df$n.male.cl <- cl2_sankey_data[match(row_df$gene, cl2_sankey_data$gene), ]$n.male.cl
row_df$n.female.cl <- cl2_sankey_data[match(row_df$gene, cl2_sankey_data$gene), ]$n.female.cl
row_df$sex <- factor(row_df$sex, levels = c("Male", "Female"))
row_df$t.cl <- factor(
  ifelse(row_df$sex == "Male", as.character(row_df$male.cl), as.character(row_df$female.cl)),
  levels = as.character(1:4)
)
row_df$n.cl <- factor(
  ifelse(row_df$sex == "Male", row_df$n.male.cl, row_df$n.female.cl),
  levels = as.character(1:4)
)
row_df$o.cl <- factor(
  ifelse(row_df$sex == "Male", as.character(row_df$female.cl), as.character(row_df$male.cl)),
  levels = as.character(1:4)
)

row_df$row.split <- factor(
  paste(row_df$sex, row_df$t.cl, row_df$n.cl),
  levels = c(
    paste("Male", rep(1:4, each = 4), rep(1:4, times = 4)),
    paste("Female", rep(1:4, each = 4), rep(1:4, times = 4))
  )
)
row_df <- row_df[order(row_df$row.split), , drop = FALSE]
rownames(row_df) <- row_df$label

make_long_matrix <- function(tissue) {
  out <- matrix(
    NA_real_,
    nrow = nrow(row_df),
    ncol = length(prediction_ids),
    dimnames = list(rownames(row_df), prediction_ids)
  )
  for (sex in c("Male", "Female")) {
    idx <- which(row_df$sex == sex)
    out[idx, ] <- fit_mats[[sex]][[tissue]][row_df$gene[idx], prediction_ids, drop = FALSE]
  }
  out
}

normal_mat <- make_long_matrix("Normal")
tumor_mat <- make_long_matrix("Tumor")
colnames(normal_mat) <- paste("Normal", prediction_ids, sep = "_")
colnames(tumor_mat) <- paste("Tumor", prediction_ids, sep = "_")

# Place tumor trajectory columns leftmost, followed by normal/reference
# trajectory columns.
plot_mat <- cbind(tumor_mat, normal_mat)
finite_plot_values <- sum(is.finite(plot_mat))
message(
  "Figure 2A AD-TMP matrix: ",
  nrow(plot_mat), " rows x ", ncol(plot_mat), " columns; ",
  finite_plot_values, " finite values; range ",
  paste(signif(range(plot_mat, na.rm = TRUE), 3), collapse = " to ")
)
if (finite_plot_values == 0) {
  stop("Figure 2A plotting matrix contains no finite AD-TMP values.", call. = FALSE)
}
col_df <- data.frame(
  id = colnames(plot_mat),
  age = rep(prediction_ages, times = 2),
  age.class = factor(rep(as.character(col_age_class), times = 2), levels = age_levels),
  tissue = factor(rep(c("Tumor", "Normal"), each = length(prediction_ids)), levels = c("Tumor", "Normal")),
  stringsAsFactors = FALSE
)

pathway_df <- data.frame(
  lapply(pathway_gene_sets[plot_pathways], function(genes) row_df$gene %in% genes),
  check.names = FALSE
)
rownames(pathway_df) <- rownames(row_df)
row_df <- data.frame(pathway_df, row_df, check.names = FALSE)

make_label_annotation <- function(row_data) {
  panel_mark_genes <- unique(mark_genes[
    mark_genes %in% row_data$gene[rowSums(row_data[, plot_pathways, drop = FALSE], na.rm = TRUE) > 0]
  ])
  mark_at <- which(row_data$gene %in% panel_mark_genes)
  mark_labels <- row_data$gene[mark_at]

  ComplexHeatmap::rowAnnotation(
    labels = ComplexHeatmap::anno_mark(
      at = mark_at,
      labels = mark_labels,
      which = "row",
      side = "right",
      labels_gp = grid::gpar(fontsize = 8),
      link_gp = grid::gpar(col = "black", lwd = 0.7),
      link_width = grid::unit(4, "mm"),
      padding = grid::unit(0.5, "mm")
    ),
    width = grid::unit(30, "mm"),
    show_annotation_name = FALSE
  )
}

make_left_annotation <- function(row_data, show_legend = TRUE) {
  ComplexHeatmap::rowAnnotation(
    df = row_data[, c("t.cl", "n.cl", "o.cl", "sex", plot_pathways), drop = FALSE],
    col = c(
      pathway_col,
      list(
        "sex" = sex_col,
        "t.cl" = cluster_col,
        "n.cl" = normal_cluster_col,
        "o.cl" = cluster_col
      )
    ),
    na_col = "gray95",
    annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    show_legend = show_legend,
    simple_anno_size = grid::unit(3, "mm")
  )
}

make_top_annotation <- function(show_legend = TRUE) {
  ComplexHeatmap::HeatmapAnnotation(
    df = col_df[, c("tissue", "age.class"), drop = FALSE],
    col = list("age.class" = age_class_col, "tissue" = tissue_col),
    annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    show_legend = show_legend,
    simple_anno_size = grid::unit(3, "mm")
  )
}

make_sex_heatmap <- function(sex, show_legends = TRUE) {
  keep <- row_df$sex == sex
  row_data <- droplevels(row_df[keep, , drop = FALSE])
  mat <- plot_mat[rownames(row_data), , drop = FALSE]

  ComplexHeatmap::Heatmap(
    mat,
    name = "Protein AD-TMP",
    col = heatmap_col,
    width = grid::unit(60, "mm"),
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    use_raster = TRUE,
    raster_quality = 5,
    na_col = "gray95",
    right_annotation = make_label_annotation(row_data),
    left_annotation = make_left_annotation(row_data, show_legend = show_legends),
    top_annotation = make_top_annotation(show_legend = show_legends),
    row_gap = grid::unit(5, "mm"),
    column_gap = grid::unit(1, "mm"),
    column_split = col_df$tissue,
    # Preserve the coordinated row ordering within each sex-specific block. The
    # combined row_df was already ordered above by sex, T-TMP group, and N-TMP
    # group.
    row_split = factor(
      paste(row_data$sex, row_data$t.cl),
      levels = paste(sex, 1:4)
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = sex,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_title = NULL,
    show_heatmap_legend = show_legends,
    heatmap_legend_param = list(
      title = "Protein AD-TMP",
      at = c(-1, -0.5, 0, 0.5, 1)
    )
  )
}

male_ht <- make_sex_heatmap("Male", show_legends = FALSE)
female_ht <- make_sex_heatmap("Female", show_legends = TRUE)

draw_long_heatmap <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 26, height = 14, useDingbats = FALSE)
  } else {
    grDevices::png(file, width = 26, height = 14, units = "in", res = 220)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 1,
    ncol = 2,
    widths = grid::unit(c(0.48, 0.52), "npc")
  )))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  ComplexHeatmap::draw(
    male_ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    newpage = FALSE,
    padding = grid::unit(c(2, 2, 2, 2), "mm")
  )
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  ComplexHeatmap::draw(
    female_ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    show_heatmap_legend = TRUE,
    show_annotation_legend = TRUE,
    newpage = FALSE,
    padding = grid::unit(c(2, 45, 2, 2), "mm")
  )
  grid::upViewport(2)
}

pdf_file <- file.path(output_dir, "Figure2A_from_temporalCPSA.pdf")
png_file <- file.path(output_dir, "Figure2A_from_temporalCPSA.png")
draw_long_heatmap(pdf_file, "pdf")
draw_long_heatmap(png_file, "png")

message("Saved ", pdf_file, " and ", png_file)
