# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

## ============================================================
## Tumor vs normal protein plots (all cancer groups)
##
## This script generates tumor-versus-normal protein abundance
## plots for selected genes across all tumor samples combined.
##
## Notes:
## - Study-derived clinical variables are harmonized to simple
##   analysis-facing names (e.g., age, Gender) for plotting
##   clarity. Original source columns are retained in the
##   clinical table.
## - Tumor protein data are read directly from the public protein
##   matrix in ../data.
## - Normal reference proteome data are loaded from temporalCPSA
##   package reference dataset.
## ============================================================

if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA before running this script.", call. = FALSE)
}

library(ggplot2)
library(patchwork)

##Helper function

mean_collapse_by_gene <- function(data_mat, gene_vec) {
  genes <- unique(gene_vec)
  
  out <- do.call(
    rbind,
    lapply(genes, function(g) {
      apply(
        data_mat[gene_vec == g, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
      )
    })
  )
  
  rownames(out) <- genes
  out
}
## ------------------------------------------------------------
## Load data
## ------------------------------------------------------------

data_dir <- "../data"

clinical.data <- temporalCPSA::ageTMP_load_clinical(data_dir)
clinical.data$id <- temporalCPSA::ageTMP_normalize_sample_ids(clinical.data$id)
if ("sample_id" %in% colnames(clinical.data)) {
  clinical.data$sample_id <- temporalCPSA::ageTMP_normalize_sample_ids(clinical.data$sample_id)
}

protein.raw <- temporalCPSA::ageTMP_load_molecular(data_dir, modality = "protein")
protein.split <- temporalCPSA::ageTMP_split_annotation_matrix(
  protein.raw,
  annotation_cols = 1:4,
  row_id = "ApprovedGeneSymbol"
)
tumor.protein.data <- mean_collapse_by_gene(
  protein.split$matrix,
  protein.split$annotation$ApprovedGeneSymbol
)

normal.reference <- temporalCPSA::ageTMP_load_normal_reference()
normal.protein.data <- normal.reference$protein$matrix
normal.protein.meta <- normal.reference$protein$sample_metadata

set.seed(123)

## ------------------------------------------------------------
## Harmonize analysis-facing clinical variable names
## ------------------------------------------------------------
## For plotting clarity, create standardized analysis aliases
## from source-specific clinical columns. 

clinical.data$age <- clinical.data$cDisc_age
clinical.data$Gender <- clinical.data$cDisc_Gender

clinical.data$cancer.group <- ifelse(
  clinical.data$cohort == "HOPE",
  as.character(clinical.data$Disc_Cancer_Group),
  "ADULT-GBM"
)

clinical.data$cancer.group <- ifelse(
  !clinical.data$Disc_Cancer_Group %in% c(
    "(DMG) Diffuse Midline Glioma",
    "ADULT-GBM",
    "(HGG) High Grade Glioma (not otherwise specified)"
  ),
  "Other",
  as.character(clinical.data$Disc_Cancer_Group)
)


clinical.data$cancer.group <- ifelse(
  clinical.data$cancer.group == "(HGG) High Grade Glioma (not otherwise specified)",
  "HGG-NOS",
  ifelse(
    clinical.data$cancer.group == "(DMG) Diffuse Midline Glioma",
    "DMG",
    ifelse(clinical.data$cancer.group == "ADULT-GBM", "ADULT-GBM", "Other")
  )
)

## ------------------------------------------------------------
## Restrict to HGG-NOS tumor samples
## ------------------------------------------------------------

nos.ids <- clinical.data[
  clinical.data$cancer.group == "HGG-NOS",
  "id"
]

## ------------------------------------------------------------
## Function: tumor vs normal plot for one gene
## ------------------------------------------------------------

get_tn_plot <- function(mygene) {
  
  ## ---- normal reference samples
  normal.df <- data.frame(
    value = as.numeric(scale(normal.protein.data[match(mygene, rownames(normal.protein.data)), ])),
    type = "normal",
    age = normal.protein.meta[match(colnames(normal.protein.data), normal.protein.meta$ID), ]$Age,
    sex = ifelse(
      normal.protein.meta[match(colnames(normal.protein.data), normal.protein.meta$ID), ]$Gender == "m",
      "Male",
      "Female"
    )
  )
  
  normal.df$age.class <- cut(
    normal.df$age,
    breaks = c(0, 15, 40),
    include.lowest = TRUE,
    labels = c("PED", "AYA")
  )
  normal.df <- normal.df[normal.df$age < 50, ]
  normal.df$value <- as.numeric(scale(normal.df$value))
  
  ## ---- tumor samples (all cancer groups combined)
  tumor.df <- data.frame(
    value = as.numeric(scale(as.numeric(tumor.protein.data[match(mygene, rownames(tumor.protein.data)), ]))),
    type = "tumor",
    age = clinical.data[match(colnames(tumor.protein.data), clinical.data$id), ]$age,
    sex = clinical.data[match(colnames(tumor.protein.data), clinical.data$id), ]$Gender
  )
  
  tumor.df <- tumor.df[tumor.df$age < 50, ]
  tumor.df$age.class <- cut(
    tumor.df$age,
    breaks = c(0, 15, 40),
    include.lowest = TRUE,
    labels = c("PED", "AYA")
  )
  tumor.df$value <- as.numeric(scale(tumor.df$value))
  
  ## ---- combine
  plot.df <- rbind(normal.df, tumor.df)
  plot.df$type <- factor(plot.df$type)
  plot.df$age.class <- factor(plot.df$age.class, levels = c("PED", "AYA"))
  plot.df$sex <- factor(plot.df$sex, levels = c("Male", "Female"))
  plot.df$label <- factor(paste(plot.df$type, plot.df$sex))
  plot.df <- plot.df[!is.na(plot.df$age.class), ]
  
  ## ---- statistical testing
  test.df <- plot.df[
    !is.na(plot.df$value) &
      !is.na(plot.df$type) &
      !is.na(plot.df$age.class) &
      !is.na(plot.df$sex),
    ,
    drop = FALSE
  ]
  pval.list <- lapply(split(test.df, list(test.df$sex, test.df$age.class), drop = TRUE), function(x) {
    if (length(unique(x$type)) <= 1) {
      return(NULL)
    }
    p <- wilcox.test(value ~ type, data = x)$p.value
    label <- if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      ""
    }
    if (label == "") {
      return(NULL)
    }
    data.frame(
      sex = x$sex[1],
      age.class = x$age.class[1],
      label = label,
      y = 2,
      stringsAsFactors = FALSE
    )
  })
  pval.df <- do.call(rbind, pval.list)
  if (is.null(pval.df)) {
    pval.df <- data.frame(
      sex = character(),
      age.class = factor(levels = levels(plot.df$age.class)),
      label = character(),
      y = numeric()
    )
  }
  pval.df$sex <- factor(pval.df$sex, levels = levels(plot.df$sex))
  pval.df$age.class <- factor(pval.df$age.class, levels = levels(plot.df$age.class))
  
  ## ---- plot
  ggplot(plot.df, aes(x = age.class, y = value, fill = label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(size = 0.5, position = position_jitterdodge(jitter.width = 0.5)) +
    geom_text(
      data = pval.df,
      aes(x = age.class, y = y, label = label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = c(
      "normal Female" = "#FFA7A7",
      "normal Male" = "#A3A3FF",
      "tumor Female" = "#CC0303",
      "tumor Male" = "#0707CF"
    )) +
    ggtitle(mygene) +
    facet_grid(cols = vars(sex)) +
    theme_bw() +
    coord_cartesian(ylim = c(-3, 3)) +
    ylab("Protein abundance") +
    xlab("Age class")
}

## ------------------------------------------------------------
## Genes to plot
## ------------------------------------------------------------

mygenes <- c("CNTN1", "MAPT", "L1CAM")
plots <- lapply(mygenes, get_tn_plot)

## ------------------------------------------------------------
## Save figure
## ------------------------------------------------------------

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
script_dir <- resolve_script_dir("Figure2D_tn_boxplot.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(output_dir, "Figure2D_protein_tn.pdf"), height = 3, width = 3 * length(mygenes))
print(wrap_plots(plots, nrow = 1, guides = "collect") & theme(legend.position = "right"))
dev.off()
