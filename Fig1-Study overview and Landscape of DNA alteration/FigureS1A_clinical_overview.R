#!/usr/bin/env Rscript

# ============================================================
# Figure S1A | Discovery Cohort Data Availability
# File: FigureS1A_clinical_overview.R
#
# Description:
#   Generates the sample-level data availability overview for the
#   Discovery cohort using the public data-availability table in
#   STable1.
#
# Inputs:
#   - ../data/STable1.xlsx, sheet Disc_DataTypeAvailabilty
#
# Outputs:
#   - FigureS1A_sample_availability.pdf
#
# Author: Weiping Ma
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

library(openxlsx)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

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
script_dir <- resolve_script_dir("FigureS1A_clinical_overview.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
# Read Excel sheet 3 (data availability table)
dat <- read.xlsx("../data/STable1.xlsx", sheet = 3)
dat <- dat[!is.na(dat$sample_id) & nzchar(as.character(dat$sample_id)), , drop = FALSE]
availability_cols <- setdiff(colnames(dat), "sample_id")
dat[availability_cols] <- lapply(dat[availability_cols], function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("TRUE", "T", "YES", "Y", "1") ~ TRUE,
    x %in% c("FALSE", "F", "NO", "N", "0") ~ FALSE,
    TRUE ~ NA
  )
})

# Quick check: column sums (all data except sample_id col)
colSums(dat[,-1])

# 2. Define sample ordering
sample_order <- dat %>%
  arrange(desc(Proteomics), 
          desc(Phosphoproteomics), 
          desc(Glycoproteomics),
          desc(WGS), 
          desc(WGS_tumor_only), 
          desc(RNAseq), 
          desc(Methylation),
          desc(snRNASeq)) %>%
  pull(sample_id)

# 3. Reshape data for plotting
dat <- melt(dat[, colnames(dat) != "snRNASeq"], 
            id.vars = "sample_id", 
            variable.name = "data_type", 
            value.name = "data_availability")

dat$data_type <- factor(dat$data_type, levels = c(
  # "snRNASeq", 
  "Methylation", "RNAseq", 
  "WGS_tumor_only", "WGS", 
  "Glycoproteomics", "Phosphoproteomics", "Proteomics"
))

dat$sample_id <- factor(dat$sample_id, levels = sample_order)

dat <- dat %>%
  mutate(label = ifelse(data_availability == TRUE, as.character(data_type), FALSE))

# 4. Define colors
color.p <- RColorBrewer::brewer.pal(n = 5, name = 'Greens')
color.r <- RColorBrewer::brewer.pal(n = 6, name = 'Blues')

# 5. Generate base plot
q <- ggplot(dat %>% filter(!data_type %in% c("WGS_tumor_only")),
            aes(sample_id, data_type, fill = label)) + 
  geom_tile(colour = "white", aes(height = 1)) + 
  ggpubr::theme_pubr() +
  scale_fill_manual(values = c(
    "FALSE" = "white", 
    "Proteomics" = color.p[5], 
    "Phosphoproteomics" = color.p[4],
    "Glycoproteomics" = color.p[3],
    "WGS" = color.r[6], 
    "RNAseq" = color.r[4], 
    "Methylation" = color.r[3]
    # ,"snRNASeq" = color.r[2]
  )) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) + 
  xlab("") + ylab("") +
  theme(legend.position = "none")

# 6. Add counts to y-axis labels
labs <- dat %>% 
  filter(label != FALSE) %>%
  group_by(label) %>% 
  summarise(n = n()) %>% 
  mutate(label2 = paste0(label, ' (n = ', n, ')')) %>%
  dplyr::select(-c(n))

# Handle WGS vs WGS_tumor_only labels
part1 <- labs %>% filter(label == "WGS") %>% pull(label2) 
part2 <- labs %>% filter(label == "WGS_tumor_only") %>% pull(label2) 

labs <- labs %>%
  filter(label != "WGS_tumor_only") %>%
  mutate(label2 = ifelse(label == "WGS", paste(part1, part2, sep = "\n"), label2)) 

# Update y-axis labels
q <- q + scale_y_discrete(
  limit = labs$label[match(setdiff(levels(dat$data_type), "WGS_tumor_only"), labs$label)],
  labels = labs$label2[match(setdiff(levels(dat$data_type), "WGS_tumor_only"), labs$label)]
)

# 7. Overlay WGS_tumor_only boxes
dat_tmp <- dat %>%
  filter(data_type == "WGS_tumor_only",
         data_availability == TRUE) %>%
  mutate(data_type = "WGS", 
         label = ifelse(label != FALSE | sample_id %in% (dat$sample_id[dat$data_type == 'WGS']), "WGS", FALSE))

q <- q + geom_tile(data = dat_tmp, aes(height = 0.5, width = 0.9), 
                   color = "black", fill = color.r[5]) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(2, 1, 1, 1, "cm")
  )

# 8. Save plot
ggsave(filename = file.path(output_dir, "FigureS1A_sample_availability.pdf"), plot = q, width = 12, height = 4)
