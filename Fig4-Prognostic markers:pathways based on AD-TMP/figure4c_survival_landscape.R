#!/usr/bin/env Rscript

# ==============================================================================
# File: Figure4C_survival_landscape.R
# ==============================================================================
#
# Title:
#   Figure 4C — Cross-cohort survival association summary heatmap
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   This script generates a heatmap summarizing sex- and age-stratified survival
#   associations for protein and RNA features in the cDiscovery/reference cohort
#   analysis. Signed local FDR-transformed survival association scores are
#   thresholded to indicate positive, negative, absent, or missing associations.
#
# Input:
#   - data/STable4.xlsx
#       Sheet: SA-Protein-cDisc-Ref
#       Sheet: SA-RNA-cDisc-Ref
#
# Output:
#   - Figure4C_survival_landscape.pdf
#
# Notes:
#   Values with absolute signed score > 1 are encoded by sign:
#     -1 = negative survival association
#      1 = positive survival association
#      0 = no association
#      3 = missing value
#
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ComplexHeatmap)
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
    return(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)))
  }
  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop("Cannot determine script directory for `", script_name, "`. Run from the script folder or with Rscript.", call. = FALSE)
}
script_dir <- resolve_script_dir("Figure4C_survival_landscape.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sex.col <- c(
  "Male" = "#0707CF",
  "Female" = "#CC0303"
)

age.class.col <- c(
  "PED" = "#B8E3B2",
  "ADO" = "#77C679",
  "YA" = "#238443",
  "ADULT" = "#00441B"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

read_survival_matrix <- function(sheet_name, feature_col = "gene") {
  
  dat <- readxl::read_excel(
    "../data/STable4.xlsx",
    sheet = sheet_name,
    na = c("NA", "")
  )
  
  female_cols <- c(
    feature_col,
    "YA.comb.signed.log10locfdr.ref.cont.female",
    "ADO.comb.signed.log10locfdr.ref.cont.female",
    "PED.comb.signed.log10locfdr.ref.cont.female"
  )
  
  male_cols <- c(
    feature_col,
    "YA.comb.signed.log10locfdr.ref.cont.male",
    "ADO.comb.signed.log10locfdr.ref.cont.male",
    "PED.comb.signed.log10locfdr.ref.cont.male"
  )
  
  fmat <- dat[, female_cols]
  mmat <- dat[, male_cols]
  
  colnames(fmat) <- gsub(
    "ref\\.|log10|\\.female",
    "",
    colnames(fmat)
  )
  
  colnames(mmat) <- gsub(
    "ref\\.|log10|\\.male",
    "",
    colnames(mmat)
  )
  
  out <- merge(
    mmat,
    fmat,
    by = feature_col,
    suffixes = c(".male", ".female"),
    all = TRUE
  )
  
  out[, -1] <- apply(
    out[, -1],
    c(1, 2),
    function(x) ifelse(is.na(x), 0, x)
  )
  
  out
}

# ------------------------------------------------------------------------------
# Read protein and RNA survival association matrices
# ------------------------------------------------------------------------------

protein.mat <- read_survival_matrix("SA-Protein-cDisc-Ref", feature_col = "gene")
rna.mat     <- read_survival_matrix("SA-RNA-cDisc-Ref", feature_col = "gene")

# ------------------------------------------------------------------------------
# Merge and threshold survival associations
# ------------------------------------------------------------------------------

mat <- merge(
  protein.mat,
  rna.mat,
  by = "gene",
  all = TRUE,
  suffixes = c(".protein", ".rna")
)

colnames(mat) <- gsub(
  "comb.signed.locfdr.cont.",
  "",
  colnames(mat)
)

rownames(mat) <- mat$gene

row.df <- data.frame(
  gene = mat$gene,
  stringsAsFactors = FALSE
)

mat <- apply(
  mat[, -1],
  c(1, 2),
  function(x) ifelse(abs(x) > 1, sign(x), 0)
)

mat <- mat[
  apply(mat, 1, function(x) sum(abs(x), na.rm = TRUE)) != 0,
]

mat <- apply(
  mat,
  c(1, 2),
  function(x) ifelse(is.na(x), 3, x)
)

# ------------------------------------------------------------------------------
# Define feature source annotations
# ------------------------------------------------------------------------------

both <- intersect(protein.mat$gene, rna.mat$gene)

protein.only <- protein.mat$gene[
  !protein.mat$gene %in% both
]

rna.only <- rna.mat$gene[
  !rna.mat$gene %in% both
]

type <- ifelse(row.df$gene %in% both, "both", NA)
type <- ifelse(row.df$gene %in% protein.only, "protein.only", type)
type <- ifelse(row.df$gene %in% rna.only, "rna.only", type)

row.df$type <- factor(
  type,
  levels = c("protein.only", "rna.only", "both")
)

row.df <- row.df[
  match(rownames(mat), row.df$gene),
]

# ------------------------------------------------------------------------------
# Define column annotations
# ------------------------------------------------------------------------------

col.df <- data.frame(
  age.class = factor(
    rep(c("YA", "ADO", "PED"), 4),
    levels = c("YA", "ADO", "PED")
  ),
  sex = factor(
    rep(rep(c("Male", "Female"), each = 3), 2),
    levels = c("Male", "Female")
  ),
  type = rep(c("protein", "rna"), each = 6)
)

age.class.col2 <- age.class.col[c("PED", "ADO", "YA")]
names(age.class.col2) <- c("PED", "ADO", "YA")

# ------------------------------------------------------------------------------
# Generate heatmap
# ------------------------------------------------------------------------------

tpsa.ht <- Heatmap(
  mat,
  name = "Surv. Assoc.",
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  row_split = factor(
    row.df$type,
    levels = c("both", "rna.only", "protein.only")
  ),
  column_split = factor(
    paste(col.df$type, col.df$sex),
    levels = c(
      "rna Female",
      "protein Female",
      "rna Male",
      "protein Male"
    )
  ),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  right_annotation = rowAnnotation(
    df = data.frame(type = row.df$type),
    col = list(
      type = c(
        "protein.only" = "#FFA500",
        "rna.only"     = "#0000FF",
        "both"         = "#A52A2A"
      )
    )
  ),
  top_annotation = HeatmapAnnotation(
    df = col.df[, c("sex", "age.class", "type")],
    col = list(
      sex = sex.col,
      type = c(
        "protein" = "#FFA500",
        "rna"     = "#0000FF",
        "both"    = "#A52A2A"
      ),
      age.class = age.class.col2
    )
  ),
  col = c(
    "-1" = "#4DAC26",
    "1"  = "#D01C8B",
    "0"  = "#F7F7F7",
    "3"  = "#E5E5E5"
  )
)

pdf(file.path(output_dir, "Figure4C_survival_landscape.pdf"), height = 8, width = 7)
draw(tpsa.ht)
dev.off()

# ------------------------------------------------------------------------------
# Summarize positive and negative associations by comparison
# ------------------------------------------------------------------------------

association.counts <- apply(
  mat,
  2,
  function(x) table(factor(x, levels = c(-1, 1)))
)

association.counts <- t(
  association.counts[, c(
    "PED.male.protein",
    "ADO.male.protein",
    "YA.male.protein",
    "PED.male.rna",
    "ADO.male.rna",
    "YA.male.rna",
    "PED.female.protein",
    "ADO.female.protein",
    "YA.female.protein",
    "PED.female.rna",
    "ADO.female.rna",
    "YA.female.rna"
  )]
)

association.counts
