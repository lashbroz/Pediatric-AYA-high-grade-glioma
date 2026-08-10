#!/usr/bin/env Rscript

# ==============================================================================
# Discovery cohort mutation cascade heatmap
# ==============================================================================
#
# Description:
#   Generates a mutation cascade heatmap for pediatric (PED) and adolescent/
#   young adult (AYA) discovery cohort high-grade glioma samples.
#
#   Mutation classes are shown by gene and sample, with genes retained if mutated
#   in >5% of discovery samples. Samples are cascade-ordered by binary mutation
#   pattern, and PED/AYA mutation-class summaries are shown as side annotations.
#
# Inputs:
#   - ../data/STable1.xlsx
#       * ClinicalTable
#       * Disc_Mutation
#
# Outputs:
#   - Figure1C_oncplot.pdf
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomic Sciences
#   Icahn School of Medicine at Mount Sinai
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
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
script_dir <- resolve_script_dir("Figure1C_oncplot.R")

data_dir <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, "..", "data"), mustWork = TRUE)
}

output_dir <- if (length(args) >= 2) args[[2]] else file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

input_workbook <- file.path(data_dir, "STable1.xlsx")

heatmap_pdf_out <- file.path(output_dir, "Figure1C_oncplot.pdf")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

standardize_mutation_class <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan")] <- NA_character_
  x <- gsub(" ", "_", x)
  x <- gsub("-", "_", x)
  x
}

mutation_to_binary <- function(x) {
  as.integer(!is.na(x) & x != "WT" & x != "0")
}

compute_binary_mutation_pct <- function(mat) {
  apply(mat, 1, function(x) {
    100 * prop.table(table(factor(x, levels = c(0, 1))))[2]
  })
}

compute_group_label <- function(mat, group) {
  apply(mat, 1, function(x) {
    pct <- tapply(x, group, function(y) {
      100 * prop.table(table(factor(y, levels = c(0, 1))))[2]
    })
    paste0(round(pct, 0), "%", collapse = " ")
  })
}

cascade_order <- function(mat) {
  if (nrow(mat) <= 1) return(mat)
  
  key <- rep("", ncol(mat))
  
  for (i in seq_len(nrow(mat) - 1)) {
    key <- paste0(key, mat[i, ])
    ord <- order(key, decreasing = TRUE)
    mat <- mat[, ord, drop = FALSE]
    key <- key[ord]
  }
  
  mat
}

make_named_palette <- function(values, palette_name = "Set3") {
  values <- sort(unique(na.omit(as.character(values))))
  if (length(values) == 0) return(character())
  
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  base <- RColorBrewer::brewer.pal(
    min(max(length(values), 3), max_colors),
    palette_name
  )
  
  cols <- colorRampPalette(base)(length(values))
  names(cols) <- values
  cols
}

mix_and_total_counts <- function(mat, idx, cats, valid_levels) {
  idx[is.na(idx)] <- FALSE
  
  if (sum(idx) == 0) {
    mix_counts <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(cats),
      dimnames = list(rownames(mat), cats)
    )
    return(list(mix_counts = mix_counts, total_pct = rep(0, nrow(mat))))
  }
  
  sub <- mat[, idx, drop = FALSE]
  
  all_counts <- t(apply(sub, 1, function(x) {
    as.numeric(table(factor(x, levels = valid_levels)))
  }))
  colnames(all_counts) <- valid_levels
  
  list(
    mix_counts = all_counts[, cats, drop = FALSE],
    total_pct = 100 * (1 - all_counts[, "WT"] / sum(idx))
  )
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

clinical_data <- as.data.frame(
  readxl::read_excel(input_workbook, sheet = "ClinicalTable")
)

mutation_status <- as.data.frame(
  readxl::read_excel(input_workbook, sheet = "Disc_Mutation")
)

# ------------------------------------------------------------------------------
# Prepare mutation and clinical matrices
# ------------------------------------------------------------------------------

mutation_genes <- mutation_status$Gene

mutation_class_matrix <- as.matrix(
  mutation_status[, setdiff(colnames(mutation_status), "Gene"), drop = FALSE]
)

rownames(mutation_class_matrix) <- mutation_genes

mutation_class_matrix <- apply(
  mutation_class_matrix,
  c(1, 2),
  standardize_mutation_class
)

mutation_class_matrix <- as.matrix(mutation_class_matrix)
rownames(mutation_class_matrix) <- mutation_genes

binary_mutation_matrix <- apply(
  mutation_class_matrix,
  c(1, 2),
  mutation_to_binary
)

binary_mutation_matrix <- as.matrix(binary_mutation_matrix)
rownames(binary_mutation_matrix) <- rownames(mutation_class_matrix)
colnames(binary_mutation_matrix) <- colnames(mutation_class_matrix)

clinical_data$Sex <- clinical_data$cDisc_Gender
clinical_data$Diagnosis <- clinical_data$cDisc_Diagnosis
clinical_data$Cancer.Group <- clinical_data$Disc_Cancer_Group
clinical_data$tumor.loc <- clinical_data$cDisc_tumor_loc
clinical_data$diagnosis.type <- clinical_data$cDisc_diagnosis_type
clinical_data$age <- as.numeric(clinical_data$cDisc_age)

rownames(clinical_data) <- clinical_data$id

col_df <- clinical_data[, c(
  "id",
  "age",
  "Diagnosis",
  "Sex",
  "Cancer.Group",
  "diagnosis.type",
  "cohort",
  "tumor.loc",
  "Disc_MutCount"
)]

rownames(col_df) <- col_df$id

col_df <- col_df[col_df$cohort == "Disc", , drop = FALSE]

shared_samples <- intersect(
  col_df$id,
  colnames(binary_mutation_matrix)
)

col_df <- col_df[shared_samples, , drop = FALSE]

binary_mutation_matrix <- binary_mutation_matrix[
  ,
  shared_samples,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  ,
  shared_samples,
  drop = FALSE
]

sex_order <- order(col_df$Sex)

binary_mutation_matrix <- binary_mutation_matrix[
  ,
  sex_order,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  ,
  sex_order,
  drop = FALSE
]

col_df <- col_df[
  colnames(binary_mutation_matrix),
  ,
  drop = FALSE
]

col_df$age.class <- factor(
  cut(
    col_df$age,
    breaks = c(0, 15, 40),
    include.lowest = TRUE,
    labels = c("PED", "AYA")
  )
)

keep_age <- !is.na(col_df$age.class)

binary_mutation_matrix <- binary_mutation_matrix[
  ,
  keep_age,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  ,
  keep_age,
  drop = FALSE
]

col_df <- col_df[
  keep_age,
  ,
  drop = FALSE
]

# ------------------------------------------------------------------------------
# Filter and order mutations
# ------------------------------------------------------------------------------

prevalence_order <- order(
  rowSums(binary_mutation_matrix, na.rm = TRUE),
  decreasing = TRUE
)

binary_mutation_matrix <- binary_mutation_matrix[
  prevalence_order,
  ,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  rownames(binary_mutation_matrix),
  ,
  drop = FALSE
]

gene_pct <- compute_binary_mutation_pct(binary_mutation_matrix)

binary_mutation_matrix <- binary_mutation_matrix[
  gene_pct > 5,
  ,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  rownames(binary_mutation_matrix),
  ,
  drop = FALSE
]

row_pct_label <- compute_group_label(
  binary_mutation_matrix,
  col_df$age.class
)

rownames(binary_mutation_matrix) <- paste(
  row_pct_label,
  rownames(binary_mutation_matrix)
)

col_df$mut.count <- suppressWarnings(
  as.numeric(col_df$Disc_MutCount)
)

computed_mut_count <- colSums(
  !is.na(mutation_class_matrix) &
    mutation_class_matrix != "WT"
)

missing_counts <- is.na(col_df$mut.count)

col_df$mut.count[missing_counts] <- computed_mut_count[
  rownames(col_df)[missing_counts]
]

binary_mutation_matrix <- cascade_order(binary_mutation_matrix)

col_df <- col_df[
  colnames(binary_mutation_matrix),
  ,
  drop = FALSE
]

mutation_class_matrix <- mutation_class_matrix[
  match(
    gsub(".* ", "", rownames(binary_mutation_matrix)),
    rownames(mutation_class_matrix)
  ),
  colnames(binary_mutation_matrix),
  drop = FALSE
]

rownames(mutation_class_matrix) <- rownames(binary_mutation_matrix)

# ------------------------------------------------------------------------------
# Define colors and annotations
# ------------------------------------------------------------------------------

mutation_colors <- c(
  "Missense_Mutation" = "#E69F00",
  "Nonsense_Mutation" = "#D55E00",
  "Frame_Shift_Del" = "#0072B2",
  "Frame_Shift_Ins" = "#56B4E9",
  "In_Frame_Del" = "#009E73",
  "In_Frame_Ins" = "#00BA87",
  "Splice_Site" = "#CC79A7",
  "Splice_Region" = "#AA4499",
  "WT" = "#EEEEEE"
)

mut_mat <- mutation_class_matrix
valid_levels <- names(mutation_colors)

mut_mat[
  !(mut_mat %in% valid_levels) | is.na(mut_mat)
] <- "WT"

age_class_col <- c(
  PED = "#7B61A8",
  AYA = "#C9A227"
)

sex_col <- c(
  "Male" = "#1F4E79",
  "Female" = "#C2185B"
)

tumor_loc_col <- c(
  "Cerebellar" = "#97005D",
  "Cortical" = "#DA8569",
  "Midline" = "#2E4B63",
  "Other/Multiple locations/NOS" = "#9EADB1"
)

diagnosis_col <- make_named_palette(col_df$Diagnosis, "Set3")
diagnosis_type_col <- make_named_palette(col_df$diagnosis.type, "Paired")
cancer_group_col <- make_named_palette(col_df$Cancer.Group, "Set2")

mut_count_max <- max(col_df$mut.count, na.rm = TRUE)

mut_count_col <- circlize::colorRamp2(
  c(0, mut_count_max),
  c("#F7FCF5", "#1B7F1B")
)


# ------------------------------------------------------------------------------
# Annotation palettes
# ------------------------------------------------------------------------------

age_class_col <- c(
  PED = "#2F92F3",
  AYA = "#21E63B"
)

sex_col <- c(
  "Male" = "#0019FF",
  "Female" = "#F00000"
)

tumor_loc_col <- c(
  "Cerebellar" = "#97005D",
  "Cortical" = "#D88B73",
  "Midline" = "#344A63",
  "Other/Multiple locations/NOS" = "#95A5AE"
)

diagnosis_type_col <- c(
  "Primary" = "#C7DD8A",
  "Progressive" = "#8B7C9E",
  "Recurrent" = "#3A2A6B",
  "Second Malignancy" = "#158C96"
)

mutation_colors <- c(
  "Missense_Mutation" = "#D89C17",
  "Nonsense_Mutation" = "#D95F02",
  "Frame_Shift_Del" = "#56B4E9",
  "Frame_Shift_Ins" = "#1F77B4",
  "In_Frame_Del" = "#1B9E77",
  "In_Frame_Ins" = "#00A65A",
  "Splice_Site" = "#CC79A7",
  "Splice_Region" = "#A55194",
  "WT" = "#E6E6E6"
)
annotation_vars <- c(
  "age.class",
  "Diagnosis",
  "diagnosis.type",
  "tumor.loc",
  "Cancer.Group",
  "Sex",
  "mut.count"
)

annotation_labels <- c(
  "age.class" = "Age class",
  "Diagnosis" = "Diagnosis",
  "diagnosis.type" = "Diagnosis type",
  "tumor.loc" = "Tumor location",
  "Cancer.Group" = "Cancer group",
  "Sex" = "Sex",
  "mut.count" = "Mutation count"
)

pct_overall <- 100 * (1 - rowSums(mut_mat == "WT") / ncol(mut_mat))
row_order <- order(pct_overall, decreasing = TRUE)
mut_mat <- mut_mat[row_order, , drop = FALSE]

age_vec <- col_df[colnames(mut_mat), "age.class"]

ped_idx <- !is.na(age_vec) & age_vec == "PED"
aya_idx <- !is.na(age_vec) & age_vec == "AYA"

cats <- setdiff(valid_levels, "WT")

ped <- mix_and_total_counts(mut_mat, ped_idx, cats, valid_levels)
aya <- mix_and_total_counts(mut_mat, aya_idx, cats, valid_levels)

ped_lab <- sprintf("%d%%", round(ped$total_pct))
aya_lab <- sprintf("%d%%", round(aya$total_pct))

max_ped <- max(rowSums(ped$mix_counts), na.rm = TRUE)
max_aya <- max(rowSums(aya$mix_counts), na.rm = TRUE)

ticks_ped <- sort(unique(c(0, pretty(c(0, max_ped), n = 4))))
ticks_aya <- sort(unique(c(0, pretty(c(0, max_aya), n = 4))))

# ------------------------------------------------------------------------------
# Build and export heatmap
# ------------------------------------------------------------------------------

top_ann <- HeatmapAnnotation(
  df = col_df[colnames(mut_mat), annotation_vars, drop = FALSE],
  col = list(
    "age.class" = age_class_col,
    "Diagnosis" = diagnosis_col,
    "diagnosis.type" = diagnosis_type_col,
    "tumor.loc" = tumor_loc_col,
    "Cancer.Group" = cancer_group_col,
    "Sex" = sex_col,
    "mut.count" = mut_count_col
  ),
  annotation_label = unname(annotation_labels[annotation_vars]),
  annotation_name_gp = gpar(fontsize = 9, fontface = "plain")
)

left_ann <- rowAnnotation(
  "PED mix (count)" = anno_barplot(
    ped$mix_counts,
    stack = TRUE,
    width = unit(32, "mm"),
    border = FALSE,
    gp = gpar(fill = mutation_colors[cats]),
    ylim = c(0, max(max_ped, 1)),
    axis_param = list(
      side = "bottom",
      at = ticks_ped,
      labels = as.character(ticks_ped),
      gp = gpar(fontsize = 7)
    )
  ),
  "PED %" = anno_text(
    ped_lab,
    gp = gpar(fontsize = 9, fontface = "bold"),
    location = 0.5,
    just = "center",
    width = unit(9, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

right_ann <- rowAnnotation(
  "AYA mix (count)" = anno_barplot(
    aya$mix_counts,
    stack = TRUE,
    width = unit(32, "mm"),
    border = FALSE,
    gp = gpar(fill = mutation_colors[cats]),
    ylim = c(0, max(max_aya, 1)),
    axis_param = list(
      side = "bottom",
      at = ticks_aya,
      labels = as.character(ticks_aya),
      gp = gpar(fontsize = 7)
    )
  ),
  "AYA %" = anno_text(
    aya_lab,
    gp = gpar(fontsize = 9, fontface = "bold"),
    location = 0.5,
    just = "center",
    width = unit(9, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

ht <- Heatmap(
  mut_mat,
  col = mutation_colors,
  na_col = "white",
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.1),
  name = "Mutation status",
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_df[colnames(mut_mat), "age.class"],
  height = nrow(mut_mat) * unit(5, "mm"),
  width = ncol(mut_mat) * unit(1, "mm"),
  top_annotation = top_ann,
  left_annotation = left_ann,
  right_annotation = right_ann
)

lgd_mix <- Legend(
  title = "Alteration type",
  labels = cats,
  legend_gp = gpar(fill = mutation_colors[cats])
)

pdf(heatmap_pdf_out, width = 20, height = 20)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_mix)
)
dev.off()

message("Wrote heatmap PDF: ", heatmap_pdf_out)
