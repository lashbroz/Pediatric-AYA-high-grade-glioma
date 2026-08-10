# ==============================================================================
# File: plot_ttmp_pathway_enrichment.R
# ==============================================================================
#
# Title:
#   Tumor Cluster Pathway Enrichment Analysis
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   This script generates a sex-stratified tumor cluster pathway enrichment
#   bar plot using pathway-level enrichment statistics from Supplementary
#   Table 2 (TTMP_Protein_Pathway sheet). Male enrichment scores are plotted
#   to the left of zero and female enrichment scores to the right.
#
# Input:
#   data/STable2.xlsx
#     └── Sheet: TTMP_Protein_Pathway
#
# Output:
#   Figure2B_barplot_tn_cluster_pathway_all.pdf
#
# Dependencies:
#   readxl
#   ggplot2
#
# Notes:
#   - Pathways are displayed in a manually curated biological order.
#   - Only pathways with FDR < 0.05 are included.
#   - Cluster colors correspond to major tumor cluster identities.
#
# ============================================================================== 
library(readxl)
library(ggplot2)

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

plotme <- read_excel(
  "../data/STable2.xlsx",
  sheet = "TTMP_Protein_Pathway"
)

# ------------------------------------------------------------------------------
# Define plotting order and colors
# ------------------------------------------------------------------------------

mylev <- c(
  "GOBP_MRNA_PROCESSING",
  "GOBP_CHROMATIN_REMODELING",
  "GOBP_NUCLEAR_TRANSPORT",
  "REACTOME_GENE_SILENCING_BY_RNA",
  "REACTOME_DEGRADATION_OF_GLI1_BY_THE_PROTEASOME",
  "REACTOME_FORMATION_OF_RNA_POL_II_ELONGATION_COMPLEX",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_NCRNA_METABOLIC_PROCESS",
  "MITO3_OXPHOS",
  "GOMF_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "REACTOME_NEURONAL_SYSTEM",
  "MITO3_MITOCHONDRIAL_CENTRAL_DOGMA",
  "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
  "GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",
  "GOBP_MITOCHONDRIAL_MEMBRANE_ORGANIZATION",
  "GOBP_REGULATION_OF_TRANSMEMBRANE_TRANSPORT",
  "GOBP_EXOCYTOSIS",
  "ADD_GLYCOPROTEIN",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "REACTOME_SELENOAMINO_ACID_METABOLISM",
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
  "REACTOME_SARS_COV_1_INFECTION",
  "ADD_XGENE"
)

cluster.col <- c(
  "1" = "#e41a1c",
  "2" = "#377eb8",
  "3" = "#3EDD51",
  "4" = "#6202ad"
)

# ------------------------------------------------------------------------------
# Format data
# ------------------------------------------------------------------------------

plotme$cl <- factor(plotme$cl)
plotme$sex <- factor(plotme$sex, levels = c("male", "female", "male.female"))

plotme$label <- paste(plotme$pathway, plotme$cl)
plotme$pure <- plotme$cl %in% c("11", "22", "33", "44")

plotme$clm <- factor(gsub("(.).", "\\1", plotme$cl))
plotme$clf <- factor(gsub(".(.)", "\\1", plotme$cl))

plotme$pathway <- factor(plotme$pathway, levels = rev(mylev))

plotme.keep <- plotme[
  plotme$fdr < 0.05 &
    plotme$sex %in% c("male", "female") &
    !is.na(plotme$pathway),
]

plotme.keep <- plotme.keep[order(plotme.keep$cl, plotme.keep$fdr), ]

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(
  plotme.keep,
  aes(
    y = pathway,
    x = ifelse(sex == "male", -signed.fdr, signed.fdr),
    fill = clm
  )
) +
  geom_bar(
    position = "dodge",
    stat = "identity"
  ) +
  scale_fill_manual(values = cluster.col, name = "Cluster") +
  geom_rect(
    xmin = -1,
    xmax = 1,
    ymin = -Inf,
    ymax = Inf,
    fill = "white",
    color = "white",
    alpha = 0.1
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    legend.position = "right"
  ) +
  xlab("Tumor cluster enrichment score")

p

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
script_dir <- resolve_script_dir("Figure2B_tn_diff_by_sex.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "Figure2B_barplot_tn_cluster_pathway_all.pdf"),
  plot = p,
  width = 10,
  height = 7,
  units = "in"
)
