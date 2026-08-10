#!/usr/bin/env Rscript

# ==============================================================================
# File: Figure4B_survival_assoc_glyco.R
# ==============================================================================
#
# Title:
#   Figure 4B — Glycopeptide survival-association heatmap
#
# Author: Nicole L. Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#
# Description:
#   Recreates Figure 4B from the public supplementary Table 4 workbook only.
#   The heatmap shows glycopeptide-level survival association scores across
#   sex, age class, and protein-adjusted/unadjusted glycopeptide analyses.
#
# Provenance:
#   - Pathway/module groups are asserted from the five rows marked
#     selected == "Yes" in STable4.xlsx, sheet SA-Glyco-Pathway-cDisc.
#   - Glycopeptide-level signed p-value colors and signed FDR significance
#     point overlays are read from STable4.xlsx, sheet SA-Glyco-Disc.
#   - Glycopeptide rows are selected programmatically using the Figure 4B
#     pathway/module gene sets plus the original glyco.int2 significance rule:
#     male PED FDR-significant in both adjusted and unadjusted models, or female
#     PED FDR-significant in either adjusted or unadjusted models. The plotted
#     row for each gene is the strongest significant membrane glycopeptide in
#     the selected pathway panel, evaluated separately by sex before union.
#   - Note: FAM234B passes this reconstructed CLUSTER_24 selection rule, although
#     it was not present in the original displayed Figure 4B panel. Table 4 does
#     not encode a final displayed-panel flag, so this script documents the
#     discrepancy rather than treating FAM234B as a failed candidate.
#
# Input:
#   data/STable4.xlsx
#     - SA-Glyco-Disc
#     - SA-Glyco-Pathway-cDisc
#
# Outputs:
#   Fig4-Prognostic markers:pathways based on AD-TMP/output/
#     - Figure4B_survival_assoc_glyco.pdf
#     - Figure4B_survival_assoc_glyco.png
#     - Figure4B_selected_glycopeptides_from_STable4.tsv
#
# ==============================================================================

required_packages <- c(
  "readxl",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_file <- gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]]))
    return(dirname(normalizePath(script_file, mustWork = TRUE)))
  }

  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }

  stop(
    "Cannot determine script directory for `", script_name,
    "`. Run from the script folder or with Rscript.",
    call. = FALSE
  )
}

script_dir <- resolve_script_dir("Figure4B_survival_assoc_glyco.R")

# Keep input/output paths aligned with the other Figure 4 table-only scripts.
input_file <- file.path(script_dir, "..", "data", "STable4.xlsx")
sheet_name <- "SA-Glyco-Disc"
output_dir <- file.path(script_dir, "output")
output_pdf <- file.path(output_dir, "Figure4B_survival_assoc_glyco.pdf")
output_png <- file.path(output_dir, "Figure4B_survival_assoc_glyco.png")
output_selected_rows <- file.path(output_dir, "Figure4B_selected_glycopeptides_from_STable4.tsv")
plot_width <- 8.2
plot_height <- 8.2
fdr_threshold <- 0.10
signed_fdr_threshold <- -log10(fdr_threshold)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

available_sheets <- readxl::excel_sheets(input_file)
required_sheets <- c(sheet_name, "SA-Glyco-Pathway-cDisc")
missing_sheets <- setdiff(required_sheets, available_sheets)
if (length(missing_sheets) > 0) {
  stop(
    "Required sheet(s) not found: ", paste(missing_sheets, collapse = ", "),
    "\nAvailable sheets: ", paste(available_sheets, collapse = ", "),
    call. = FALSE
  )
}

raw_data <- readxl::read_excel(input_file, sheet = sheet_name, col_types = "text")
pathway_data <- readxl::read_excel(input_file, sheet = "SA-Glyco-Pathway-cDisc", col_types = "text")

required_columns <- c(
  "glycoppetide",
  "membrane.loc",
  "PED.comb.signed.p.female.adj",
  "ADO.comb.signed.p.female.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.pmale.adj",
  "ADO.comb.signed.pmale.adj",
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.p.female.noadj",
  "ADO.comb.signed.p.female.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj",
  "PED.comb.signed.pmale.noadj",
  "ADO.comb.signed.pmale.noadj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj"
)

missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0) {
  stop(
    "Required column(s) missing from ", sheet_name, ": ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

numeric_columns <- setdiff(required_columns, c("glycoppetide", "membrane.loc"))
for (column in numeric_columns) {
  raw_data[[column]] <- suppressWarnings(as.numeric(raw_data[[column]]))
}

raw_data$gene <- sub("_.*", "", raw_data$glycoppetide)
raw_data$glycan <- sub(".*-", "", raw_data$glycoppetide)
raw_data$display_label <- paste(raw_data$gene, raw_data$glycan, sep = " - ")

pathway_name_map <- c(
  "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT" =
    "Regulation of\nNervous System\nDevelopment",
  "CLUSTER_24" = "TMP Group MF-24",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION" =
    "Extracellular\nMatrix\nOrganization",
  "CLUSTER_13" = "TMP Group MF-13",
  "HALLMARK_COAGULATION" = "Coagulation"
)

# These five selected pathway/module rows are the provenance for the Figure 4B
# row groups. This is a guardrail: if Table 4 changes, the script should fail
# loudly instead of silently drawing a different biological panel.
selected_pathways <- pathway_data$pathway[pathway_data$selected == "Yes"]
expected_pathways <- names(pathway_name_map)
if (!setequal(selected_pathways, expected_pathways)) {
  stop(
    "The selected pathways in SA-Glyco-Pathway-cDisc do not match the Figure 4B pathway set.\n",
    "Expected: ", paste(expected_pathways, collapse = ", "), "\n",
    "Observed selected == Yes: ", paste(selected_pathways, collapse = ", "),
    call. = FALSE
  )
}

# The pathway sheet stores pathway-level statistics but not member genes.
# These gene sets are the Figure 4B pathway/module definitions from the original
# pathway analysis source, restricted here to glycopeptide-covered genes needed
# to derive the displayed panel. Glycopeptide rows below are selected by criteria,
# not by a hand-written list of final marker labels.
pathway_gene_sets <- list(
  "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT" = c(
    "BCAN", "EPHB3", "ITGB1", "L1CAM", "MEGF8", "NPTN", "NTRK2",
    "PLXNA4", "PLXNB2", "PLXNB3", "PLXNC1", "PTPRZ1", "THY1"
  ),
  "CLUSTER_24" = c(
    # FAM234B is retained in the module definition because it passes the
    # reconstructed criteria, but it was not part of the original displayed
    # Figure 4B panel; see provenance note above.
    "CNTN1", "FAM234B", "OMG"
  ),
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION" = c(
    "ADAM17", "BCAN", "CD44", "CD47", "COL6A2", "ITGA1", "ITGA3",
    "ITGAV", "ITGB1", "ITGB3", "LAMC1", "NCAM1", "NID2", "PLOD1",
    "PLOD3", "PXDN"
  ),
  "CLUSTER_13" = c(
    "CD276", "LAMP1"
  ),
  "HALLMARK_COAGULATION" = c(
    "LAMP2", "LRP1"
  )
)

missing_pathway_sets <- setdiff(expected_pathways, names(pathway_gene_sets))
if (length(missing_pathway_sets) > 0) {
  stop(
    "Internal pathway gene-set definition(s) missing: ",
    paste(missing_pathway_sets, collapse = ", "),
    call. = FALSE
  )
}

# Candidate glycopeptides are deduplicated in the original pathway selection
# order, so genes present in multiple selected pathways keep the first pathway
# assignment used in the Figure 4B source workflow.
pathway_assignment_order <- c(
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT",
  "CLUSTER_24",
  "CLUSTER_13",
  "HALLMARK_COAGULATION"
)

male_ped_adjusted <- raw_data$glycoppetide[
  abs(raw_data[["PED.comb.signed.fdrmale.adj"]]) > signed_fdr_threshold
]
male_ped_unadjusted <- raw_data$glycoppetide[
  abs(raw_data[["PED.comb.signed.fdrmale.noadj"]]) > signed_fdr_threshold
]
female_ped_adjusted <- raw_data$glycoppetide[
  abs(raw_data[["PED.comb.signed.fdr.female.adj"]]) > signed_fdr_threshold
]
female_ped_unadjusted <- raw_data$glycoppetide[
  abs(raw_data[["PED.comb.signed.fdr.female.noadj"]]) > signed_fdr_threshold
]

glyco_int2 <- unique(c(
  intersect(male_ped_adjusted, male_ped_unadjusted),
  female_ped_adjusted,
  female_ped_unadjusted
))

candidate_rows <- do.call(
  rbind,
  lapply(pathway_assignment_order, function(pathway_id) {
    candidate_index <- which(
      raw_data$gene %in% pathway_gene_sets[[pathway_id]] &
        raw_data$glycoppetide %in% glyco_int2 &
        raw_data$membrane.loc == "Yes"
    )

    data.frame(
      source_pathway = pathway_id,
      row_index = candidate_index,
      glycopeptide = raw_data$glycoppetide[candidate_index],
      gene = raw_data$gene[candidate_index],
      display_label = raw_data$display_label[candidate_index],
      stringsAsFactors = FALSE
    )
  })
)

candidate_rows <- candidate_rows[
  match(unique(candidate_rows$glycopeptide), candidate_rows$glycopeptide),
  ,
  drop = FALSE
]

if (nrow(candidate_rows) == 0) {
  stop("No glycopeptide rows passed the Figure 4B selection criteria.", call. = FALSE)
}

fdr_columns <- c(
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj"
)

p_columns <- c(
  "PED.comb.signed.pmale.adj",
  "ADO.comb.signed.pmale.adj",
  "PED.comb.signed.p.female.adj",
  "ADO.comb.signed.p.female.adj",
  "PED.comb.signed.pmale.noadj",
  "ADO.comb.signed.pmale.noadj",
  "PED.comb.signed.p.female.noadj",
  "ADO.comb.signed.p.female.noadj"
)

male_fdr_columns <- c(
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj"
)

female_fdr_columns <- c(
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj"
)

select_best_glycopeptide_by_gene <- function(fdr_column_names) {
  fdr_score <- apply(
    abs(as.data.frame(raw_data[candidate_rows$row_index, fdr_column_names, drop = FALSE])),
    1,
    function(x) max(x, na.rm = TRUE)
  )
  fdr_score[is.infinite(fdr_score)] <- NA_real_

  sex_rows <- candidate_rows
  sex_rows$max_abs_signedlog10fdr <- fdr_score
  sex_rows <- sex_rows[
    !is.na(sex_rows$max_abs_signedlog10fdr) &
      sex_rows$max_abs_signedlog10fdr > signed_fdr_threshold,
    ,
    drop = FALSE
  ]
  sex_rows <- sex_rows[order(abs(sex_rows$max_abs_signedlog10fdr), decreasing = TRUE), , drop = FALSE]

  sex_rows[match(unique(sex_rows$gene), sex_rows$gene), , drop = FALSE]
}

male_best_rows <- select_best_glycopeptide_by_gene(male_fdr_columns)
female_best_rows <- select_best_glycopeptide_by_gene(female_fdr_columns)
best_glycopeptides <- unique(c(male_best_rows$glycopeptide, female_best_rows$glycopeptide))

selected_candidates <- candidate_rows[candidate_rows$glycopeptide %in% best_glycopeptides, , drop = FALSE]
selected_candidates$source_pathway <- factor(selected_candidates$source_pathway, levels = names(pathway_name_map))
selected_candidates <- selected_candidates[
  order(selected_candidates$source_pathway, selected_candidates$gene, selected_candidates$display_label),
  ,
  drop = FALSE
]

selected_index <- selected_candidates$row_index
plot_rows <- raw_data[selected_index, , drop = FALSE]
plot_rows$source_pathway <- as.character(selected_candidates$source_pathway)
plot_rows$group <- factor(
  unname(pathway_name_map[plot_rows$source_pathway]),
  levels = unname(pathway_name_map)
)
plot_rows$published_label <- plot_rows$display_label

max_fdr_score <- apply(
  abs(as.data.frame(plot_rows[, fdr_columns, drop = FALSE])),
  1,
  function(x) max(x, na.rm = TRUE)
)
max_fdr_score[is.infinite(max_fdr_score)] <- NA_real_

plot_rows$selected_by_male <- plot_rows$glycoppetide %in% male_best_rows$glycopeptide
plot_rows$selected_by_female <- plot_rows$glycoppetide %in% female_best_rows$glycopeptide
plot_rows$max_abs_signedlog10fdr <- max_fdr_score

# Audit table: documents the criterion-selected glycopeptides and their pathway
# membership source, including the sex-specific rule that selected each row.
utils::write.table(
  data.frame(
    selected_order = seq_len(nrow(plot_rows)),
    source_pathway = plot_rows$source_pathway,
    pathway_group = gsub("\n", " ", as.character(plot_rows$group)),
    display_label = plot_rows$published_label,
    table4_glycopeptide = plot_rows$glycoppetide,
    gene = plot_rows$gene,
    glycan = plot_rows$glycan,
    selected_by_male = plot_rows$selected_by_male,
    selected_by_female = plot_rows$selected_by_female,
    max_abs_signedlog10fdr = plot_rows$max_abs_signedlog10fdr,
    membrane.loc = plot_rows$membrane.loc,
    stringsAsFactors = FALSE
  ),
  file = output_selected_rows,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Tile colors show signed -log10(FDR) values. Black point overlays mark
# entries where the corresponding signed -log10(FDR) reaches FDR < 0.10.
value_columns <- c(
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj"
)

significance_columns <- c(
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj"
)

heatmap_matrix <- as.matrix(plot_rows[, value_columns, drop = FALSE])
significance_matrix <- as.matrix(plot_rows[, significance_columns, drop = FALSE])
significance_matrix <- ifelse(
  abs(significance_matrix) >= signed_fdr_threshold,
  significance_matrix,
  NA_real_
)

column_labels <- c("PED", "ADO", "PED", "ADO", "PED", "ADO", "PED", "ADO")
colnames(heatmap_matrix) <- column_labels
colnames(significance_matrix) <- column_labels
rownames(heatmap_matrix) <- plot_rows$published_label
rownames(significance_matrix) <- plot_rows$published_label

column_metadata <- data.frame(
  Adjustment = factor(
    c(rep("Adj.", 4), rep("No Adj.", 4)),
    levels = c("Adj.", "No Adj.")
  ),
  Sex = factor(
    rep(c("M", "F"), each = 2, times = 2),
    levels = c("M", "F")
  ),
  `Age class` = factor(
    rep(c("PED", "ADO"), times = 4),
    levels = c("PED", "ADO")
  ),
  check.names = FALSE
)

column_split <- factor(
  paste(column_metadata$Adjustment, column_metadata$Sex),
  levels = unique(paste(column_metadata$Adjustment, column_metadata$Sex))
)

plot_rows$sialylated <- grepl("S[1-9]", plot_rows$glycan)
plot_rows$fucosylated <- grepl("F[1-9]", plot_rows$glycan)
plot_rows$high_mannose <- grepl("^N2H[5-9]F0S0G0$", plot_rows$glycan)

pathway_colors <- c(
  "Regulation of\nNervous System\nDevelopment" = "#b7d900",
  "TMP Group MF-24" = "#c018d9",
  "Extracellular\nMatrix\nOrganization" = "#2ecc71",
  "TMP Group MF-13" = "#2146d0",
  "Coagulation" = "#d7191c"
)

sex_colors <- c("M" = "#0000cc", "F" = "#d00000")
age_class_colors <- c("PED" = "#bfe7bf", "ADO" = "#6fbd6a")
adjustment_colors <- c("Adj." = "#f7ad0b", "No Adj." = "#8e4d07")
absent_glycan_color <- "#eeeeee"
glycan_symbol_colors <- c(
  Sial = "#ff1a1a",
  Fuco = "#b238ff",
  HM = "#13e95d"
)

heatmap_colors <- circlize::colorRamp2(
  c(-2, 0, 2),
  rev(colorRampPalette(RColorBrewer::brewer.pal(5, "PiYG"))(3))
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = column_metadata,
  col = list(
    Adjustment = adjustment_colors,
    Sex = sex_colors,
    `Age class` = age_class_colors
  ),
  annotation_name_side = "left",
  show_legend = TRUE
)

left_annotation <- ComplexHeatmap::rowAnnotation(
  Pathway = plot_rows$group,
  Sial = ComplexHeatmap::anno_simple(
    rep("absent", nrow(plot_rows)),
    which = "row",
    col = c(absent = absent_glycan_color),
    pch = ifelse(plot_rows$sialylated, 23, 16),
    pt_size = grid::unit(2.8, "mm"),
    pt_gp = grid::gpar(
      fill = ifelse(plot_rows$sialylated, glycan_symbol_colors[["Sial"]], absent_glycan_color),
      col = ifelse(plot_rows$sialylated, "black", absent_glycan_color),
      lwd = 0.8
    )
  ),
  Fuco = ComplexHeatmap::anno_simple(
    rep("absent", nrow(plot_rows)),
    which = "row",
    col = c(absent = absent_glycan_color),
    pch = ifelse(plot_rows$fucosylated, 24, 16),
    pt_size = grid::unit(2.8, "mm"),
    pt_gp = grid::gpar(
      fill = ifelse(plot_rows$fucosylated, glycan_symbol_colors[["Fuco"]], absent_glycan_color),
      col = ifelse(plot_rows$fucosylated, "black", absent_glycan_color),
      lwd = 0.8
    )
  ),
  HM = ComplexHeatmap::anno_simple(
    rep("absent", nrow(plot_rows)),
    which = "row",
    col = c(absent = absent_glycan_color),
    pch = ifelse(plot_rows$high_mannose, 21, 16),
    pt_size = grid::unit(2.8, "mm"),
    pt_gp = grid::gpar(
      fill = ifelse(plot_rows$high_mannose, glycan_symbol_colors[["HM"]], absent_glycan_color),
      col = ifelse(plot_rows$high_mannose, "black", absent_glycan_color),
      lwd = 0.8
    )
  ),
  col = list(
    Pathway = pathway_colors
  ),
  annotation_name_rot = 90,
  annotation_name_gp = grid::gpar(fontsize = 7),
  simple_anno_size = grid::unit(3.4, "mm"),
  show_legend = TRUE
)

glyco_heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Glyco-peptide\nSurv. Assoc.",
  col = heatmap_colors,
  na_col = "grey88",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = plot_rows$group,
  column_split = column_split,
  row_gap = grid::unit(c(1.5, 1.5, 1.5, 1.5), "mm"),
  column_gap = grid::unit(c(1.2, 2.6, 1.2), "mm"),
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.6),
  top_annotation = top_annotation,
  left_annotation = left_annotation,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 6.9),
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 7),
  column_names_rot = 90,
  column_title = NULL,
  width = grid::unit(ncol(heatmap_matrix) * 4.8, "mm"),
  height = grid::unit(nrow(heatmap_matrix) * 4.1, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(significance_matrix[i, j])) {
      grid::grid.points(
        x,
        y,
        pch = 16,
        size = grid::unit(1.05, "mm"),
        gp = grid::gpar(col = "black")
      )
    }
  }
)

fdr_point_legend <- ComplexHeatmap::Legend(
  title = "Significance",
  labels = paste0("FDR < ", fdr_threshold),
  type = "points",
  pch = 16,
  size = grid::unit(1.05, "mm"),
  legend_gp = grid::gpar(col = "black"),
  title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 7)
)

glycan_type_legend <- ComplexHeatmap::Legend(
  title = "Glycan Type",
  labels = c("Sialylated", "Fucosylated", "High Mannose"),
  type = "points",
  pch = c(23, 24, 21),
  size = grid::unit(3.2, "mm"),
  legend_gp = grid::gpar(
    fill = unname(glycan_symbol_colors[c("Sial", "Fuco", "HM")]),
    col = "black"
  ),
  title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 7)
)

message(
  "Using signed FDR columns for both tile color and FDR < ",
  fdr_threshold, " point overlays."
)
message(
  "Figure row groups are the five rows marked selected == Yes in ",
  "SA-Glyco-Pathway-cDisc. Glycopeptide rows are selected from pathway/module ",
  "gene sets using the original glyco.int2 PED FDR rule and strongest ",
  "significant membrane glycopeptide per gene; all plotted association values ",
  "are read from ", input_file, "."
)

draw_panel <- function() {
  grid::grid.newpage()
  grid::grid.text(
    "Glyco-peptide Survival Assoc.",
    x = grid::unit(0.57, "npc"),
    y = grid::unit(0.975, "npc"),
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
  ComplexHeatmap::draw(
    glyco_heatmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(fdr_point_legend, glycan_type_legend),
    newpage = FALSE,
    padding = grid::unit(c(10, 4, 4, 4), "mm")
  )
}

grDevices::cairo_pdf(output_pdf, width = plot_width, height = plot_height)
draw_panel()
grDevices::dev.off()

grDevices::png(output_png, width = plot_width, height = plot_height, units = "in", res = 300, type = "cairo")
draw_panel()
grDevices::dev.off()

message("Saved ", output_pdf)
message("Saved ", output_png)
message("Saved ", output_selected_rows)
