#!/usr/bin/env Rscript

# ==============================================================================
# File: Figure6K_sex_biased_glycopeptide_pathway_enrichment.R
# ==============================================================================
#
# Title:
#   Figure 6K — Sex-biased glycopeptide pathway enrichment analysis
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   This standalone script reproduces the Figure 6K glycopeptide pathway
#   enrichment bar plot directly from the published supplementary tables.
#   The script is a minimal refactor of the original Figure 6K plotting
#   logic from `code_dump/sial_score_rev.R`, replacing upstream global
#   objects with explicit supplementary-table inputs and local definitions.
#
# Biological context:
#   The analysis highlights sex-biased glycopeptide pathway enrichment
#   associated with C2/M2-F2 tumor biology.
#
# Input:
#   ../data/STable6.xlsx
#     └── Sheet: C2-Specific-Sex-Bias-Pathway
#
# Output:
#   Figure6K_glycopeptide_pathway_enrichment.pdf
#
# Dependencies:
#   ggplot2
#   readxl
#
# Notes:
#   - Reproduces the original manuscript pathway filtering logic.
#   - Retains the curated Figure 6K pathway set and ordering.
#   - Uses signed log(FDR) values for directional enrichment plotting.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
})

input_file <- file.path("../data", "STable6.xlsx")
input_sheet <- "C2-Specific-Sex-Bias-Pathway"
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
script_dir <- resolve_script_dir("Figure6K_sex_biased_glycopeptide_pathway_enrichment.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "Figure6K_glycopeptide_pathway_enrichment.pdf")

# Curated pathway set used in the original Figure 6K plotting block. The order
# below is the manuscript curation order; the plotted y-axis order follows the
# original code by sorting these selected rows by signed FDR.
selected_pathways <- c(
  "REACTOME_COLLAGEN_FORMATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "GOMF_GLYCOSYLTRANSFERASE_ACTIVITY",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_MTORC1_SIGNALING"
)

# Original palette was colorspace::darken("blue", .2) and
# colorspace::darken("red", .2). Hex colors are specified here to avoid loading
# the full historical helper stack.
sex_col <- c("Male" = "#0000CC", "Female" = "#CC0000")

format_pathway_label <- function(x) {
  unlist(lapply(x, function(y) {
    y <- unlist(strsplit(as.character(y), "_"))
    y <- paste(y[1], tools::toTitleCase(tolower(paste(y[-1], collapse = " "))))
    y <- gsub(" Notch1 ", " NOTCH1 ", y)
    y <- gsub(" Mrna ", " mRNA ", y)
    y <- gsub(" Atp ", " ATP ", y)
    y <- gsub(" Adora2b ", " ADORA2B ", y)
    y <- gsub(" Ampa ", " AMPA ", y)
    y <- gsub(" Trna ", " tRNA ", y)
    y <- gsub(" Upr", " UPR", y)
    y <- gsub(" Dna ", " DNA ", y)
    y <- gsub(" And ", " and ", y)
    y <- gsub(" In ", " in ", y)
    y <- gsub(" Of ", " of ", y)
    y <- gsub(" Tca ", " TCA ", y)
    y <- gsub(" Egfr", " EGFR ", y)
    y <- gsub(" Rna ", " RNA ", y)
    y <- gsub(" Mhc ", " MHC ", y)
    y <- gsub(" Rna$", " RNA", y)
    y <- gsub(" Rrna", " rRNA", y)
    y <- gsub(" Tdp43 ", " TDP-43 ", y)
    y <- gsub(" Snrna ", " snRNA ", y)
    y <- gsub(" Npc ", " NPC ", y)
    y <- gsub(" Srp ", " SRP ", y)
    y <- gsub(" Of ", " of ", y)
    y <- gsub(" Ii ", " II ", y)
    y <- gsub(" By ", " by ", y)
    y <- gsub(" To ", " to ", y)
    y <- gsub(" And ", " and ", y)
    y <- gsub(" Nmda ", " NMDA ", y)
    y <- gsub(" Gtpase ", " GTPase ", y)
    y <- gsub(" Rhobtb2 ", " rhoBTB2 ", y)
    y <- gsub(" Gpcr", " GPCR", y)
    y
  }))
}

pathway_data <- as.data.frame(read_excel(input_file, sheet = input_sheet))

required_columns <- c("pathway", "logsignfdr", "type")
missing_columns <- setdiff(required_columns, colnames(pathway_data))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s) in ", input_file, " / ", input_sheet, ": ",
    paste(missing_columns, collapse = ", ")
  )
}

pathway_data$logsignfdr <- as.numeric(pathway_data$logsignfdr)

# Minimal refactor of the original plotting logic in sial_score_rev.R:
#   1. order all pathway rows by signed FDR,
#   2. retain glycopeptide/glycoprotein rows,
#   3. retain significant rows with abs(logsignfdr) > 1,
#   4. retain the curated Figure 6K pathways,
#   5. plot signed FDR values with sex-specific coloring.
plot_data <- pathway_data[order(pathway_data$logsignfdr, decreasing = FALSE), ]
plot_data$pathway <- factor(plot_data$pathway, levels = unique(plot_data$pathway))
plot_data$pathway.label <- factor(
  format_pathway_label(plot_data$pathway),
  levels = format_pathway_label(unique(plot_data$pathway))
)
plot_data <- plot_data[grepl("glyco", plot_data$type), ]
plot_data <- plot_data[abs(plot_data$logsignfdr) > 1, ]
plot_data <- plot_data[plot_data$pathway %in% selected_pathways, ]
plot_data <- plot_data[order(plot_data$logsignfdr), ]
plot_data$pathway.label <- factor(plot_data$pathway.label, levels = unique(plot_data$pathway.label))
plot_data$sex <- ifelse(plot_data$logsignfdr > 1, "Male", "Female")

if (nrow(plot_data) == 0) {
  stop("No Figure 6K rows remained after applying the original filters.")
}

figure6k <- ggplot(plot_data, aes(y = pathway.label, x = logsignfdr, fill = sex)) +
  geom_rect(
    data = plot_data[1, ],
    xmin = -1,
    xmax = 1,
    ymin = -Inf,
    ymax = Inf,
    col = "gray30",
    fill = "white",
    lty = 2,
    alpha = 0
  ) +
  geom_bar(stat = "identity", color = "white", alpha = 1) +
  theme_bw() +
  scale_fill_manual(values = sex_col) +
  ylab("")

ggsave(
  filename = output_file,
  plot = figure6k,
  width = 5.5,
  height = 2,
  units = "in",
  useDingbats = FALSE
)

message("Saved Figure 6K to: ", output_file)
