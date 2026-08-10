################################################################################
# Figure: Figure 4D
# Script name: Figure4D_survival_assoc_heatmap.R
# Purpose: Generate the Figure 4D protein survival-association heatmap
#          from supplementary table inputs.
#
# Inputs:
#   data/STable4.xlsx
#     - SA-Protein-cDisc-Ref
#     - SA-Protein-Pathway-cDisc-Ref
#
# Outputs:
#   Figure4D_protein_survival_overlap.pdf
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
################################################################################

set.seed(1)

required_packages <- c(
  "readxl", "ComplexHeatmap", "circlize",
  "RColorBrewer", "colorspace", "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(colorspace)
library(grid)

input_file <- file.path("../data", "STable4.xlsx")
protein_sheet <- "SA-Protein-cDisc-Ref"
pathway_sheet <- "SA-Protein-Pathway-cDisc-Ref"

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
script_dir <- resolve_script_dir("Figure4D_survival_assoc_heatmap.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "Figure4D_protein_survival_overlap.pdf")

selected_pathways <- c(
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT",
  "MITO3_MITOCHONDRIAL_RIBOSOME",
  "CLUSTER_13",
  "CLUSTER_23",
  "GOBP_RNA_SPLICING",
  "GOMF_CHROMATIN_BINDING",
  "GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY",
  "CLUSTER_14",
  "GOMF_ENDOPEPTIDASE_ACTIVITY",
  "GOBP_EXOCYTOSIS",
  "GOBP_RIBOSE_PHOSPHATE_METABOLIC_PROCESS",
  "GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H",
  "GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS",
  "REACTOME_NEURONAL_SYSTEM",
  "CLUSTER_44",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "HALLMARK_COAGULATION",
  "REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE",
  "GOMF_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "MITO3_OXPHOS",
  "CLUSTER_24"
)

figure4d_protein_pathway_map <- data.frame(
  protein = c(
    "GRN", "LRCH4",
    "NGLY1",
    "SRSF1",
    "MORC2", "NUDT21",
    "ING4", "LPXN",
    "BCAM",
    "CLPP",
    "RAB12",
    "ACACA",
    "CBR3",
    "PHKG2",
    "AP2S1", "NCALD",
    "CTSA", "STAT6"
  ),
  pathway = c(
    "CLUSTER_13", "CLUSTER_13",
    "CLUSTER_23",
    "GOBP_RNA_SPLICING",
    "GOMF_CHROMATIN_BINDING", "GOMF_CHROMATIN_BINDING",
    "GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY", "GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY",
    "CLUSTER_14",
    "GOMF_ENDOPEPTIDASE_ACTIVITY",
    "GOBP_EXOCYTOSIS",
    "GOBP_RIBOSE_PHOSPHATE_METABOLIC_PROCESS",
    "GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H",
    "GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS",
    "REACTOME_NEURONAL_SYSTEM", "REACTOME_NEURONAL_SYSTEM",
    "CLUSTER_44", "CLUSTER_44"
  ),
  stringsAsFactors = FALSE
)

selected_proteins <- figure4d_protein_pathway_map$protein

check_required_columns <- function(df, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(
      table_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

as_numeric_matrix <- function(df) {
  out <- as.data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))))
  as.matrix(out)
}

is_yes <- function(x) {
  toupper(trimws(as.character(x))) == "YES"
}

get_pathway_plot_label <- function(x) {
  vapply(x, function(y) {
    pieces <- unlist(strsplit(as.character(y), "_"))
    label <- paste(
      pieces[1],
      paste(tools::toTitleCase(tolower(pieces[-1])), collapse = " ")
    )
    
    replacements <- c(
      " Notch1 " = " NOTCH1 ",
      " Mrna " = " mRNA ",
      " Atp " = " ATP ",
      " Adora2b " = " ADORA2B ",
      " Ampa " = " AMPA ",
      " Trna " = " tRNA ",
      " Upr" = " UPR",
      " Dna " = " DNA ",
      " And " = " and ",
      " In " = " in ",
      " Of " = " of ",
      " Tca " = " TCA ",
      " Egfr" = " EGFR ",
      " Rna " = " RNA ",
      " Mhc " = " MHC ",
      " Rna$" = " RNA",
      " Rrna" = " rRNA",
      " Tdp43 " = " TDP-43 ",
      " Snrna " = " snRNA ",
      " Npc " = " NPC ",
      " Srp " = " SRP ",
      " Ii " = " II ",
      " By " = " by ",
      " To " = " to ",
      " Nmda " = " NMDA ",
      " Gtpase " = " GTPase ",
      " Rhobtb2 " = " rhoBTB2 ",
      " Gpcr" = " GPCR"
    )
    
    for (pattern in names(replacements)) {
      label <- gsub(pattern, replacements[[pattern]], label)
    }
    
    label
  }, character(1))
}

make_survival_matrix <- function(df, row_column, columns) {
  check_required_columns(df, c(row_column, columns), row_column)
  mat <- as_numeric_matrix(df[, columns, drop = FALSE])
  rownames(mat) <- as.character(df[[row_column]])
  colnames(mat) <- c("PED", "ADO", "YA", "PED", "ADO", "YA")
  mat
}

make_pathway_ref_only_mask <- function(ref_mat, cdisc_mat) {
  ref_sig <- abs(ref_mat) > 1
  cdisc_sig <- abs(cdisc_mat) > 1
  ref_sig & !cdisc_sig
}

make_pathway_concordance_mask <- function(ref_mat, cdisc_mat) {
  ref_sig <- abs(ref_mat) > 1
  cdisc_sig <- abs(cdisc_mat) > 1
  ref_sig & cdisc_sig & sign(ref_mat) == sign(cdisc_mat)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

protein_data <- readxl::read_excel(input_file, sheet = protein_sheet, na = c("", "NA"))
pathway_data <- readxl::read_excel(input_file, sheet = pathway_sheet, na = c("", "NA"))

protein_required_columns <- c(
  "gene",
  "protein.concordant",
  "PED.comb.signed.log10p.ref.male",
  "ADO.comb.signed.log10p.ref.male",
  "YA.comb.signed.log10p.ref.male",
  "PED.comb.signed.log10locfdr.ref.cont.male",
  "ADO.comb.signed.log10locfdr.ref.cont.male",
  "YA.comb.signed.log10locfdr.ref.cont.male",
  "PED.comb.signed.log10p.cdisc.male",
  "ADO.comb.signed.log10p.cdisc.male",
  "YA.comb.signed.log10p.cdisc.male",
  "PED.comb.signed.log10p.ref.female",
  "ADO.comb.signed.log10p.ref.female",
  "YA.comb.signed.log10p.ref.female",
  "PED.comb.signed.log10locfdr.ref.cont.female",
  "ADO.comb.signed.log10locfdr.ref.cont.female",
  "YA.comb.signed.log10locfdr.ref.cont.female",
  "PED.comb.signed.log10p.cdisc.female",
  "ADO.comb.signed.log10p.cdisc.female",
  "YA.comb.signed.log10p.cdisc.female"
)

pathway_required_columns <- c(
  "pathway",
  "protein.concordant.pathway",
  "selected",
  "PED.signedlog10fdr.male.ref",
  "ADO.signedlog10fdr.male.ref",
  "YA.signedlog10fdr.male.ref",
  "PED.signedlog10fdr.female.ref",
  "ADO.signedlog10fdr.female.ref",
  "YA.signedlog10fdr.female.ref",
  "PED.signedlog10fdr.male.cdisc",
  "ADO.signedlog10fdr.male.cdisc",
  "YA.signedlog10fdr.male.cdisc",
  "PED.signedlog10fdr.female.cdisc",
  "ADO.signedlog10fdr.female.cdisc",
  "YA.signedlog10fdr.female.cdisc"
)

check_required_columns(protein_data, protein_required_columns, protein_sheet)
check_required_columns(pathway_data, pathway_required_columns, pathway_sheet)

if (anyDuplicated(protein_data$gene) > 0) {
  stop("Protein sheet contains duplicated gene identifiers.", call. = FALSE)
}

if (anyDuplicated(pathway_data$pathway) > 0) {
  stop("Protein-pathway sheet contains duplicated pathway identifiers.", call. = FALSE)
}

selected_from_sheet <- as.character(pathway_data$pathway[is_yes(pathway_data$selected)])
missing_selected_pathways <- setdiff(selected_pathways, selected_from_sheet)

if (length(missing_selected_pathways) > 0) {
  stop(
    "The following hard-coded Figure 4D pathways are not marked selected in ",
    pathway_sheet, ": ",
    paste(missing_selected_pathways, collapse = ", "),
    call. = FALSE
  )
}

missing_proteins <- setdiff(selected_proteins, protein_data$gene)

if (length(missing_proteins) > 0) {
  stop(
    "The following Figure 4D proteins were not found in the protein matrix: ",
    paste(missing_proteins, collapse = ", "),
    call. = FALSE
  )
}

pathway_ref_columns <- c(
  "PED.signedlog10fdr.male.ref",
  "ADO.signedlog10fdr.male.ref",
  "YA.signedlog10fdr.male.ref",
  "PED.signedlog10fdr.female.ref",
  "ADO.signedlog10fdr.female.ref",
  "YA.signedlog10fdr.female.ref"
)

pathway_cdisc_columns <- c(
  "PED.signedlog10fdr.male.cdisc",
  "ADO.signedlog10fdr.male.cdisc",
  "YA.signedlog10fdr.male.cdisc",
  "PED.signedlog10fdr.female.cdisc",
  "ADO.signedlog10fdr.female.cdisc",
  "YA.signedlog10fdr.female.cdisc"
)

pathway_ref_mat <- make_survival_matrix(pathway_data, "pathway", pathway_ref_columns)
pathway_cdisc_mat <- make_survival_matrix(pathway_data, "pathway", pathway_cdisc_columns)

pathway_ref_mat <- pathway_ref_mat[selected_pathways, , drop = FALSE]
pathway_cdisc_mat <- pathway_cdisc_mat[selected_pathways, , drop = FALSE]

pathway_ref_only <- make_pathway_ref_only_mask(pathway_ref_mat, pathway_cdisc_mat)
pathway_concordant <- make_pathway_concordance_mask(pathway_ref_mat, pathway_cdisc_mat)

protein_ref_p_columns <- c(
  "PED.comb.signed.log10p.ref.male",
  "ADO.comb.signed.log10p.ref.male",
  "YA.comb.signed.log10p.ref.male",
  "PED.comb.signed.log10p.ref.female",
  "ADO.comb.signed.log10p.ref.female",
  "YA.comb.signed.log10p.ref.female"
)

protein_ref_locfdr_columns <- c(
  "PED.comb.signed.log10locfdr.ref.cont.male",
  "ADO.comb.signed.log10locfdr.ref.cont.male",
  "YA.comb.signed.log10locfdr.ref.cont.male",
  "PED.comb.signed.log10locfdr.ref.cont.female",
  "ADO.comb.signed.log10locfdr.ref.cont.female",
  "YA.comb.signed.log10locfdr.ref.cont.female"
)

protein_cdisc_columns <- c(
  "PED.comb.signed.log10p.cdisc.male",
  "ADO.comb.signed.log10p.cdisc.male",
  "YA.comb.signed.log10p.cdisc.male",
  "PED.comb.signed.log10p.cdisc.female",
  "ADO.comb.signed.log10p.cdisc.female",
  "YA.comb.signed.log10p.cdisc.female"
)

protein_ref_p_mat <- make_survival_matrix(protein_data, "gene", protein_ref_p_columns)
protein_ref_locfdr_mat <- make_survival_matrix(protein_data, "gene", protein_ref_locfdr_columns)
protein_cdisc_mat <- make_survival_matrix(protein_data, "gene", protein_cdisc_columns)

protein_ref_p_mat <- protein_ref_p_mat[selected_proteins, , drop = FALSE]
protein_ref_locfdr_mat <- protein_ref_locfdr_mat[selected_proteins, , drop = FALSE]
protein_cdisc_mat <- protein_cdisc_mat[selected_proteins, , drop = FALSE]

protein_ref_symbol_mat <- apply(
  protein_ref_locfdr_mat,
  c(1, 2),
  function(x) ifelse(abs(x) > 1, x, NA)
)

protein_cdisc_symbol_mat <- apply(
  protein_cdisc_mat,
  c(1, 2),
  function(x) ifelse(abs(x) > 2, x, NA)
)

rownames(protein_ref_symbol_mat) <- rownames(protein_ref_locfdr_mat)
colnames(protein_ref_symbol_mat) <- colnames(protein_ref_locfdr_mat)

rownames(protein_cdisc_symbol_mat) <- rownames(protein_cdisc_mat)
colnames(protein_cdisc_symbol_mat) <- colnames(protein_cdisc_mat)

protein_ref_only <- !is.na(protein_ref_symbol_mat) & is.na(protein_cdisc_symbol_mat)
protein_concordant <- !is.na(protein_ref_symbol_mat) & !is.na(protein_cdisc_symbol_mat)

column_data <- data.frame(
  sex = factor(c(rep("Male", 3), rep("Female", 3)), levels = c("Male", "Female")),
  age.class = factor(c("PED", "ADO", "YA", "PED", "ADO", "YA"), levels = c("PED", "ADO", "YA"))
)

sex_col <- c(
  "Female" = colorspace::darken("red", 0.2),
  "Male" = colorspace::darken("blue", 0.2)
)

age_class_col <- c(
  "PED" = "#D9F0A3",
  "ADO" = "#78C679",
  "YA" = "#238443"
)

pathway_col <- colorspace::darken(rainbow(length(selected_pathways)), 0.2)
names(pathway_col) <- selected_pathways

top_annotation <- HeatmapAnnotation(
  df = column_data,
  col = list(sex = sex_col, age.class = age_class_col)
)

pathway_color_fun <- circlize::colorRamp2(
  c(-5, 0, 5),
  rev(colorRampPalette(RColorBrewer::brewer.pal(5, "PiYG"))(3))
)

protein_color_fun <- circlize::colorRamp2(
  c(-3, 0, 3),
  rev(colorRampPalette(RColorBrewer::brewer.pal(5, "PiYG"))(3))
)

pathway_plot_mat <- pathway_ref_mat
rownames(pathway_plot_mat) <- get_pathway_plot_label(rownames(pathway_plot_mat))

protein_plot_mat <- protein_ref_p_mat

protein_row_split <- factor(
  figure4d_protein_pathway_map$pathway,
  levels = selected_pathways
)

protein_row_annotation <- rowAnnotation(
  pathway = factor(
    figure4d_protein_pathway_map$pathway,
    levels = selected_pathways
  ),
  col = list(pathway = pathway_col),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

pathway_ht <- Heatmap(
  pathway_plot_mat,
  name = "pathway\nRef FDR",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = pathway_color_fun,
  column_split = column_data$sex,
  row_gap = unit(0, "mm"),
  show_row_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  border = TRUE,
  height = nrow(pathway_plot_mat) * unit(4.5, "mm"),
  width = ncol(pathway_plot_mat) * unit(4.5, "mm"),
  heatmap_legend_param = list(direction = "horizontal"),
  left_annotation = rowAnnotation(
    pathway = factor(selected_pathways, levels = selected_pathways),
    col = list(pathway = pathway_col),
    show_annotation_name = FALSE,
    show_legend = FALSE
  ),
  top_annotation = top_annotation,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(pathway_ref_only[i, j]) && pathway_ref_only[i, j]) {
      grid.points(
        x, y,
        pch = 16,
        size = unit(1.5, "mm"),
        gp = gpar(fill = "black", col = "black")
      )
    }
    
    if (!is.na(pathway_concordant[i, j]) && pathway_concordant[i, j]) {
      grid.points(
        x, y,
        pch = 15,
        size = unit(3, "mm"),
        gp = gpar(fill = "black", col = "black")
      )
    }
  }
)

protein_ht <- Heatmap(
  protein_plot_mat,
  name = "protein\nRef p",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_split = protein_row_split,
  row_title = NULL,
  row_gap = unit(0.8, "mm"),
  col = protein_color_fun,
  column_split = column_data$sex,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  border = TRUE,
  height = nrow(protein_plot_mat) * unit(4.5, "mm"),
  width = ncol(protein_plot_mat) * unit(4.5, "mm"),
  heatmap_legend_param = list(direction = "horizontal"),
  left_annotation = protein_row_annotation,
  top_annotation = top_annotation,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(protein_ref_only[i, j]) && protein_ref_only[i, j]) {
      grid.points(
        x, y,
        pch = 21,
        size = unit(1.5, "mm"),
        gp = gpar(fill = "black", col = "black")
      )
    }
    
    if (!is.na(protein_concordant[i, j]) && protein_concordant[i, j]) {
      grid.points(
        x, y,
        pch = 22,
        size = unit(3, "mm"),
        gp = gpar(fill = "black", col = "black")
      )
    }
  }
)

pdf(output_file, width = 12, height = 8.5)

pushViewport(viewport(layout = grid.layout(
  nrow = 1,
  ncol = 2,
  widths = unit(c(0.58, 0.42), "npc")
)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(
  pathway_ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  newpage = FALSE
)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(
  protein_ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  newpage = FALSE
)
popViewport()

popViewport()
dev.off()

message("Saved Figure 4D protein survival-association heatmap to: ", output_file)
