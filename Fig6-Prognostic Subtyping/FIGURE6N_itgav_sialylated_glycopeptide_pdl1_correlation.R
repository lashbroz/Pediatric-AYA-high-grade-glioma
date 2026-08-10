#!/usr/bin/env Rscript

# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

# -----------------------------------------------------------------------------
# Figure 6N: ITGAV sialylated glycopeptide vs PD-L1/CD274 protein abundance
#
# Purpose:
#   Reproduce the Figure 6N scatter/correlation panel from repository raw data.
#
# Original datac provenance:
#   In the analysis dump, datac was built from subtype.df, clinical/survival
#   metadata, cDisc proteome values, and cDisc glycopeptide v2 values. The panel
#   then used age < 40 tumors with complete CD274 protein and
#   ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0 glycopeptide measurements. This script
#   does not read/import datac; it rebuilds the plotted columns from the raw
#   repository tables using those same variable definitions.
#
# Output:
#   output/Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.pdf
#   output/Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.png
#   output/Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation_stats.tsv
# -----------------------------------------------------------------------------

required_packages <- c("ggplot2", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

library(ggplot2)
library(readxl)

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- sub("^--file=", "", script_arg[[1]])
    script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }

  stop(
    "Cannot determine script directory for `",
    script_name,
    "`. Run from the script folder or with Rscript.",
    call. = FALSE
  )
}

z_score <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

extract_feature_row <- function(data, feature_column, feature_id, sample_ids) {
  feature_index <- match(feature_id, data[[feature_column]])
  if (is.na(feature_index)) {
    stop("Feature not found: ", feature_id, call. = FALSE)
  }

  available_ids <- intersect(sample_ids, colnames(data))
  if (length(available_ids) == 0) {
    stop("No requested sample IDs were found for feature: ", feature_id, call. = FALSE)
  }

  values <- as.numeric(data[feature_index, available_ids, drop = TRUE])
  names(values) <- available_ids
  values
}

format_r <- function(r_value) {
  out <- formatC(round(r_value, 2), format = "f", digits = 2)
  sub("0$", "", out)
}

format_p <- function(p_value) {
  if (p_value < 0.001) {
    return("p < 0.001")
  }
  paste0("p = ", formatC(p_value, format = "f", digits = 3))
}

cor_summary <- function(data, group_label) {
  test <- cor.test(data$itgav.sial.v2, data$cd274, method = "pearson")
  data.frame(
    group = group_label,
    n = nrow(data),
    r = unname(test$estimate),
    p = test$p.value,
    label = paste0("R = ", format_r(unname(test$estimate)), ", ", format_p(test$p.value)),
    stringsAsFactors = FALSE
  )
}

script_dir <- resolve_script_dir("FIGURE6N_itgav_sialylated_glycopeptide_pdl1_correlation.R")
repo_root <- dirname(script_dir)
data_dir <- file.path(repo_root, "data")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

input_stable6 <- file.path(data_dir, "STable6.xlsx")
input_clinical <- file.path(data_dir, "STable1.xlsx")
input_proteome <- file.path(data_dir, "cDisc_proteome_imputed_data_09152023.tsv")
input_glyco <- file.path(data_dir, "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv")

output_pdf <- file.path(output_dir, "Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.pdf")
output_png <- file.path(output_dir, "Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.png")
output_stats <- file.path(output_dir, "Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation_stats.tsv")

for (input_file in c(input_stable6, input_clinical, input_proteome, input_glyco)) {
  if (!file.exists(input_file)) {
    stop("Missing required input file: ", input_file, call. = FALSE)
  }
}

target_glycopeptide <- "ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0"
target_protein <- "CD274"

subtypes <- as.data.frame(readxl::read_excel(input_stable6, sheet = "Subtype-cDisc"))
subtypes <- subtypes[, c("id", "protein.subtype")]
colnames(subtypes) <- c("id", "subtype")

clinical <- as.data.frame(
  readxl::read_excel(input_clinical, sheet = "ClinicalTable"),
  check.names = FALSE
)
clinical <- clinical[, c("id", "cDisc_age", "cDisc_Gender")]
colnames(clinical) <- c("id", "age", "sex_clinical")

plot_data <- merge(subtypes, clinical, by = "id", all.x = TRUE, sort = FALSE)
plot_data$sex <- ifelse(grepl("^M", plot_data$subtype), "Male", "Female")
plot_data$sex <- ifelse(is.na(plot_data$sex), plot_data$sex_clinical, plot_data$sex)

proteome <- read.delim(
  input_proteome,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

glyco <- read.delim(
  input_glyco,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_ids <- plot_data$id

cd274_values <- extract_feature_row(
  proteome,
  "ApprovedGeneSymbol",
  target_protein,
  sample_ids
)

itgav_values <- extract_feature_row(
  glyco,
  "Gene.Sequence",
  target_glycopeptide,
  sample_ids
)

common_ids <- intersect(names(cd274_values), names(itgav_values))
measurement_data <- data.frame(
  id = common_ids,
  cd274_raw = cd274_values[common_ids],
  itgav_sial_raw = itgav_values[common_ids],
  stringsAsFactors = FALSE
)

# These z-scores mirror the datac column definitions before the age filter.
measurement_data$cd274 <- z_score(measurement_data$cd274_raw)
measurement_data$itgav.sial.v2 <- z_score(measurement_data$itgav_sial_raw)

plot_data <- merge(plot_data, measurement_data, by = "id", all.x = FALSE, sort = FALSE)
plot_data <- plot_data[
  plot_data$age < 40 &
    is.finite(plot_data$cd274) &
    is.finite(plot_data$itgav.sial.v2) &
    !is.na(plot_data$sex),
]
plot_data$sex <- factor(plot_data$sex, levels = c("Female", "Male"))

if (nrow(plot_data) < 3) {
  stop("Too few complete samples to draw Figure 6N.", call. = FALSE)
}

stats <- rbind(
  cor_summary(plot_data, "Both"),
  cor_summary(plot_data[plot_data$sex == "Female", ], "Female"),
  cor_summary(plot_data[plot_data$sex == "Male", ], "Male")
)
write.table(stats, output_stats, sep = "\t", row.names = FALSE, quote = FALSE)

sex_col <- c("Female" = "#CC0000", "Male" = "#0000CC", "Both" = "black")
sex_shape <- c("Female" = 8, "Male" = 18, "Both" = NA)
line_type <- c("Female" = "solid", "Male" = "solid", "Both" = "22")

label_df <- data.frame(
  label = stats$label,
  x = c(0.15, -0.45, -0.25),
  y = c(2.55, 2.20, 1.86),
  color = c("black", sex_col[["Female"]], sex_col[["Male"]]),
  stringsAsFactors = FALSE
)

p <- ggplot(plot_data, aes(x = itgav.sial.v2, y = cd274)) +
  geom_point(aes(color = sex, shape = sex), size = 2.2, stroke = 0.8) +
  geom_smooth(
    aes(color = sex, linetype = sex),
    method = "lm",
    se = FALSE,
    linewidth = 0.9
  ) +
  geom_smooth(
    aes(color = "Both", linetype = "Both", group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 0.9
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7,
    fontface = "italic",
    color = label_df$color
  ) +
  scale_color_manual(
    name = "Sex",
    values = sex_col,
    breaks = c("Female", "Male", "Both")
  ) +
  scale_shape_manual(
    name = "Sex",
    values = sex_shape,
    breaks = c("Female", "Male", "Both")
  ) +
  scale_linetype_manual(
    name = "Sex",
    values = line_type,
    breaks = c("Female", "Male", "Both")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(8, 18, NA),
        linetype = c("solid", "solid", "22"),
        linewidth = c(0.9, 0.9, 0.9)
      )
    ),
    shape = "none",
    linetype = "none"
  ) +
  labs(
    title = "PD-L1 Protein vs ITGAV Glyco-peptide",
    subtitle = "ANTTQPGIVEGGQVLK - N3H5F1S1G0",
    x = "ITGAV Glycopeptide (N3H5F1S1G0)\nAbundance",
    y = "CD274/PD-L1 Protein Abundance"
  ) +
  coord_cartesian(xlim = c(-2.6, 2.7), ylim = c(-2.2, 2.8), clip = "off") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5, color = "gray20"),
    legend.position = "bottom",
    legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.7),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray88", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray93", linewidth = 0.35),
    plot.margin = margin(10, 10, 10, 10)
  )

pdf(output_pdf, width = 4.25, height = 5.8, useDingbats = FALSE)
print(p)
dev.off()

png(output_png, width = 1275, height = 1740, res = 300)
print(p)
dev.off()

message("Wrote Figure 6N PDF: ", output_pdf)
message("Wrote Figure 6N PNG: ", output_png)
message("Wrote Figure 6N stats: ", output_stats)
message(
  "Samples plotted: ",
  nrow(plot_data),
  " (Female = ",
  sum(plot_data$sex == "Female"),
  ", Male = ",
  sum(plot_data$sex == "Male"),
  ")"
)
message("Correlations:")
for (i in seq_len(nrow(stats))) {
  message(
    "  ",
    stats$group[[i]],
    ": n=",
    stats$n[[i]],
    ", r=",
    formatC(stats$r[[i]], format = "f", digits = 3),
    ", p=",
    formatC(stats$p[[i]], format = "g", digits = 4)
  )
}
