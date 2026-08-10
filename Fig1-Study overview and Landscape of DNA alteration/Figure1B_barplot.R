#!/usr/bin/env Rscript

# ==============================================================================
# File: Figure1B_barplot.R
# ==============================================================================
#
# Title:
#   Figure 1B — Subtype distribution across developmental age groups
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   Generates a stacked bar plot showing tumor subtype proportions in pediatric
#   (PED) and adolescent/young adult (AYA) tumors using the public clinical
#   table in STable1.
#
# Input file:
#   - data/STable1.xlsx, sheet ClinicalTable
#
# Output file:
#   - Figure1B_barplot_subtype.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------------------
# Input / output
# ------------------------------------------------------------------------------

clinical_file <- "../data/STable1.xlsx"
clinical_sheet <- "ClinicalTable"
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
script_dir <- resolve_script_dir("Figure1B_barplot.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "Figure1B_barplot_subtype.pdf")

if (!file.exists(clinical_file)) {
  stop("Clinical workbook not found: ", clinical_file)
}
if (!clinical_sheet %in% readxl::excel_sheets(clinical_file)) {
  stop("Clinical sheet not found in STable1.xlsx: ", clinical_sheet)
}

# ------------------------------------------------------------------------------
# Read clinical data
# ------------------------------------------------------------------------------

clinical_table <- readxl::read_excel(
  clinical_file,
  sheet = clinical_sheet
) %>%
  as.data.frame()

# ------------------------------------------------------------------------------
# Preprocess
# ------------------------------------------------------------------------------

clinical_table <- clinical_table %>%
  mutate(
    age_class = ifelse(
      cDisc_age_class_name_derived %in% c("ADO", "YA"),
      "AYA",
      as.character(cDisc_age_class_name_derived)
    ),
    cancer_group = ifelse(
      flag == 1,
      "Flag",
      as.character(Disc_Cancer_Group)
    ),
    cancer_group = ifelse(
      cancer_group %in% c("NA", "", "NaN"),
      NA,
      cancer_group
    )
  )

plot_df <- clinical_table %>%
  filter(
    age_class %in% c("PED", "AYA"),
    !is.na(cancer_group)
  )

# ------------------------------------------------------------------------------
# Compute proportions
# ------------------------------------------------------------------------------

plot_data <- plot_df %>%
  count(age_class, cancer_group, name = "n") %>%
  group_by(age_class) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup() %>%
  complete(
    age_class = c("PED", "AYA"),
    cancer_group,
    fill = list(n = 0, prop = 0)
  ) %>%
  mutate(
    age_class = factor(age_class, levels = c("PED", "AYA")),
    label = ifelse(prop > 0, paste0(round(prop * 100), "%"), "")
  )

# ------------------------------------------------------------------------------
# Colors
# ------------------------------------------------------------------------------

color_palette <- c(
  "(DHG) Diffuse Hemispheric Glioma" = "#A5772B",
  "(DMG) Diffuse Midline Glioma" = "#1B3360",
  "(HGG) High Grade Glioma (not otherwise specified)" = "#4393C3",
  "(IHG) Infantile Hemispheric Glioma" = "#E62A8A",
  "(PXA) Pleomorphic Xanthoastrocytoma" = "#660821",
  "Flag" = "gray80"
)

missing_colors <- setdiff(unique(plot_data$cancer_group), names(color_palette))

if (length(missing_colors) > 0) {
  stop(
    "Missing colors for subtype(s): ",
    paste(missing_colors, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = age_class,
    y = prop,
    fill = cancer_group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    linewidth = 0.1
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = color_palette) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  labs(
    x = "Age class",
    y = "Proportion",
    fill = "Subtype"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

ggsave(
  filename = output_file,
  plot = p,
  width = 11,
  height = 3,
  useDingbats = FALSE
)

message("Wrote: ", output_file)
