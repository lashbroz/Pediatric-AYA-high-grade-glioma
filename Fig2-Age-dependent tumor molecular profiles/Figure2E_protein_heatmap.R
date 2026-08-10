# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

###############################################################
# HOPE Pediatric/AYA High-Grade Glioma Study
# Pathway-Level Tumor-Normal and Sex-Bias Heatmaps
#
# This script generates pathway-level heatmaps based on
# precomputed enrichment statistics provided in Supplementary
# Tables associated with the manuscript.
#
# Input:
#   data/STable2.xlsx
#
# Required sheets:
#   - Diff_Pathway_Protein
#   - SexBias_Pathway_Protein
#
# Output:
#   - dist_pathway_Protein TMT_tn.tsv
#   - sex_bias_pathway_Protein TMT_tn.tsv
#   - difference_heatmap_Protein TMT.pdf
#
# Notes:
#   - Heatmap fill corresponds to signed -log10(FDR)
#   - Black dots indicate nominal significance
#     (abs(signed -log10(P)) > 1)
#   - This script reproduces visualization panels from the
#     HOPE pediatric/AYA HGG proteogenomic analysis.
#
# Nicole L. Tignor, PhD
# Icahn School of Medicine at Mount Sinai
###############################################################

## -----------------------------------------------------------
## Load required libraries
## -----------------------------------------------------------

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

get_green <- function(n) {
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-(1:2)])(n)
}

## -----------------------------------------------------------
## Analysis settings
## -----------------------------------------------------------

multi <- "Protein TMT"

sex.col <- c(
  Male = "#0707CF",
  Female = "#CC0303"
)

## Selected pathways displayed in the heatmaps
mypathways <- c(
  "GOMF_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "MITO3_OXPHOS",
  "GOBP_NEURON_DEVELOPMENT",
  "GOBP_CELL_CELL_SIGNALING",
  "REACTOME_NEURONAL_SYSTEM",
  "GOMF_CALMODULIN_BINDING",
  "24"
)

## -----------------------------------------------------------
## Read pathway enrichment results from supplementary tables
## -----------------------------------------------------------

dist.data <- readxl::read_xlsx(
  "../data/STable2.xlsx",
  sheet = "Diff_Pathway_Protein",
  na = "NA"
)

sex_bias.data <- readxl::read_xlsx(
  "../data/STable2.xlsx",
  sheet = "SexBias_Pathway_Protein",
  na = "NA"
)

dist.data <- as.data.frame(dist.data)
sex_bias.data <- as.data.frame(sex_bias.data)

## -----------------------------------------------------------
## Keep required columns
## -----------------------------------------------------------

dist.data <- dist.data[, c(
  "Pathway",
  "ALL.Male.SignedLog10P",
  "ALL.Female.SignedLog10P",
  "ALL.Male.SignedLog10FDR",
  "ALL.Female.SignedLog10FDR"
)]

sex_bias.data <- sex_bias.data[, c(
  "Pathway",
  "SexBias_PED_Normal_SignedLog10P",
  "SexBias_ADO_Normal_SignedLog10P",
  "SexBias_YA_Normal_SignedLog10P",
  "SexBias_PED_Tumor_SignedLog10P",
  "SexBias_ADO_Tumor_SignedLog10P",
  "SexBias_YA_Tumor_SignedLog10P",
  "SexBias_PED_Normal_SignedLog10FDR",
  "SexBias_ADO_Normal_SignedLog10FDR",
  "SexBias_YA_Normal_SignedLog10FDR",
  "SexBias_PED_Tumor_SignedLog10FDR",
  "SexBias_ADO_Tumor_SignedLog10FDR",
  "SexBias_YA_Tumor_SignedLog10FDR"
)]

## -----------------------------------------------------------
## Tumor-normal pathway difference heatmap
## -----------------------------------------------------------

dist.data.plot <- dist.data[
  match(mypathways, dist.data$Pathway),
  ,
  drop = FALSE
]

rownames(dist.data.plot) <- dist.data.plot$Pathway

## Heatmap fill matrix
data1 <- as.matrix(dist.data.plot[, c(
  "ALL.Male.SignedLog10FDR",
  "ALL.Female.SignedLog10FDR"
)])

## Dot overlay matrix
data2 <- as.matrix(dist.data.plot[, c(
  "ALL.Male.SignedLog10P",
  "ALL.Female.SignedLog10P"
)])

## Display dots only for nominally significant pathways
data2 <- apply(
  data2,
  c(1, 2),
  function(x) ifelse(abs(x) > 1, x, NA)
)

rownames(data1) <- mypathways
rownames(data2) <- mypathways

colnames(data1) <- c("Male", "Female")
colnames(data2) <- c("Male", "Female")

## Column annotations
col.df <- data.frame(
  sex = factor(
    c("Male", "Female"),
    levels = c("Male", "Female")
  )
)

rownames(col.df) <- colnames(data1)

## Generate heatmap
p.dist0 <- Heatmap(
  as.matrix(data1),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  col = circlize::colorRamp2(
    c(0, 5),
    c("white", "#B8860B")
  ),
  
  show_column_names = FALSE,
  border = TRUE,
  
  column_split = col.df$sex,
  column_title = NULL,
  
  heatmap_legend_param = list(
    direction = "horizontal"
  ),
  
  name = "TN Euclidean Difference\n(Signed -log10 FDR)",
  
  height = nrow(data1) * unit(5, "mm"),
  width  = ncol(data1) * unit(5, "mm"),
  
  top_annotation = HeatmapAnnotation(
    df = col.df[, "sex", drop = FALSE],
    col = list(sex = sex.col),
    annotation_name_side = "left"
  ),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    if (!is.na(data2[i, j])) {
      
      grid.points(
        x,
        y,
        pch = 16,
        size = unit(2, "mm"),
        gp = gpar(
          fill = "black",
          col = "black"
        )
      )
      
    }
    
  }
  
)

## -----------------------------------------------------------
## Sex-bias pathway heatmap
## -----------------------------------------------------------

p.cols <- c(
  "SexBias_PED_Normal_SignedLog10P",
  "SexBias_ADO_Normal_SignedLog10P",
  "SexBias_YA_Normal_SignedLog10P",
  "SexBias_PED_Tumor_SignedLog10P",
  "SexBias_ADO_Tumor_SignedLog10P",
  "SexBias_YA_Tumor_SignedLog10P"
)

fdr.cols <- c(
  "SexBias_PED_Normal_SignedLog10FDR",
  "SexBias_ADO_Normal_SignedLog10FDR",
  "SexBias_YA_Normal_SignedLog10FDR",
  "SexBias_PED_Tumor_SignedLog10FDR",
  "SexBias_ADO_Tumor_SignedLog10FDR",
  "SexBias_YA_Tumor_SignedLog10FDR"
)

sex_bias.plot <- sex_bias.data[
  match(mypathways, sex_bias.data$Pathway),
  ,
  drop = FALSE
]

rownames(sex_bias.plot) <- sex_bias.plot$Pathway

## Heatmap fill matrix
data1 <- as.matrix(sex_bias.plot[, p.cols])

## Dot overlay matrix
data2 <- as.matrix(sex_bias.plot[, fdr.cols])

data2 <- apply(
  data2,
  c(1, 2),
  function(x) ifelse(abs(x) > 1, x, NA)
)

rownames(data1) <- mypathways
rownames(data2) <- mypathways

## Column annotations
col.df <- data.frame(
  tissue = c(rep("Normal", 3), rep("Tumor", 3)),
  age.class = rep(c("PED", "ADO", "YA"), 2)
)

col.df$tissue <- factor(
  col.df$tissue,
  levels = c("Normal", "Tumor")
)

col.df$age.class <- factor(
  col.df$age.class,
  levels = c("PED", "ADO", "YA")
)

col.df$split <- col.df$tissue

## Age-class colors
age.class.col <- get_green(4)[1:3]
names(age.class.col) <- c("PED", "ADO", "YA")

## Generate heatmap
pb <- Heatmap(
  
  as.matrix(data1),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  col = circlize::colorRamp2(
    c(-5, 0, 5),
    c(
      sex.col[["Female"]],
      "white",
      sex.col[["Male"]]
    )
  ),
  
  show_column_names = FALSE,
  border = TRUE,
  
  column_split = col.df$split,
  
  heatmap_legend_param = list(
    direction = "horizontal"
  ),
  
  height = nrow(data1) * unit(5, "mm"),
  width = ncol(data1) * unit(5, "mm"),
  
  name = "MF Difference\n(Signed -log10 P)",
  
  top_annotation = HeatmapAnnotation(
    df = col.df[, c("age.class", "tissue")],
    col = list(
      tissue = c(
        "Normal" = "purple",
        "Tumor" = "orange"
      ),
      age.class = age.class.col
    )
  ),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    if (!is.na(data2[i, j])) {
      
      grid.points(
        x,
        y,
        pch = 16,
        size = unit(2, "mm"),
        gp = gpar(
          fill = "black",
          col = "black"
        )
      )
      
    }
    
  }
  
)

## -----------------------------------------------------------
## Export combined heatmap figure
## -----------------------------------------------------------

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
script_dir <- resolve_script_dir("Figure2E_protein_heatmap.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pdf(
  file.path(output_dir, "Figure2E_Protein_Pathway_Heatmaps.pdf"),
  width = 10,
  height = 5
)

draw(
  p.dist0 + pb,
  heatmap_legend_side = "bottom"
)

dev.off()
