#!/usr/bin/env Rscript

# ============================================================
# Mutation Survival Associations in the External Reference Cohort
# File: mutation_survival_external_reference.R
#
# Description:
#   Evaluates mutation-survival associations in the external
#   reference cohort using multivariable Cox proportional hazards
#   models across age-defined subsets.
#
#   Models include recurrently observed mutation variables together
#   with age, sex, and diagnosis covariates. Results are summarized
#   as signed -log10(Bonferroni-adjusted p-values) and visualized
#   in a bubble plot for PED and AYA groups.
#
# Input:
#   - data/STable1.xlsx, sheet: Ref_ClinicalTable
#
# Outputs:
#   - mutation_survival_association_table.tsv
#   - Figure1E_reference_survival_assoc.pdf
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(survival)
  library(ggplot2)
  library(scales)
})

# ---- I/O ----

input_file <- "../data/STable1.xlsx"
input_sheet <- "Ref_ClinicalTable"

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
script_dir <- resolve_script_dir("Figure1E_reference_assoc.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_table <- file.path(output_dir, "mutation_survival_association_table.tsv")
output_plot <- file.path(output_dir, "Figure1E_reference_survival_assoc.pdf")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ---- Load reference clinical and mutation data ----

ref.data <- data.frame(
  read_xlsx(
    input_file,
    sheet = input_sheet,
    na = c("NA", "", "NaN")
  ),
  check.names = FALSE
)

if (nrow(ref.data) == 0) {
  stop("Reference cohort data not found or empty.")
}

required_cols <- c(
  "id",
  "Gender",
  "os.status",
  "days",
  "age",
  "age_class_name",
  "First.Diagnosis",
  "Diag"
)

missing_cols <- setdiff(required_cols, colnames(ref.data))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

ref.data$id <- as.character(ref.data$id)
row.names(ref.data) <- ref.data$id
ref.data$Sex <- ref.data$Gender

# ---- Extract mutation matrix ----

mut_cols <- grep("_mut$", colnames(ref.data), value = TRUE)

if (length(mut_cols) == 0) {
  stop("No mutation columns ending in '_mut' were found.")
}

vali.mut.data <- ref.data[, mut_cols, drop = FALSE]
row.names(vali.mut.data) <- ref.data$id

ref.meta <- ref.data[
  ,
  setdiff(colnames(ref.data), mut_cols),
  drop = FALSE
]

ref.data <- data.frame(
  ref.meta,
  vali.mut.data[match(ref.meta$id, row.names(vali.mut.data)), , drop = FALSE],
  check.names = FALSE
)

# ---- Basic filtering ----

ref.data <- ref.data[
  !is.na(ref.data$os.status) &
    !is.na(ref.data$days) &
    !is.na(ref.data$age),
  ,
  drop = FALSE
]

ref.data <- ref.data[
  ref.data$age <= 80,
  ,
  drop = FALSE
]

ref.data$age_class_name <- factor(
  ref.data$age_class_name,
  levels = c("PED", "ADO", "YA", "ADULT", "SEN")
)

# ---- Restrict to samples with observed mutation data ----

vali.mut.data <- vali.mut.data[
  row.names(vali.mut.data) %in% ref.data$id,
  ,
  drop = FALSE
]

vali.mut.cdata <- vali.mut.data[
  apply(!is.na(vali.mut.data), 1, sum) > 2,
  ,
  drop = FALSE
]

vali.mut.cdata <- vali.mut.cdata[
  ,
  apply(!is.na(vali.mut.cdata), 2, sum) == nrow(vali.mut.cdata),
  drop = FALSE
]

vali.cdata <- ref.data[
  ref.data$id %in% row.names(vali.mut.cdata),
  ,
  drop = FALSE
]

vali.cdata <- data.frame(
  vali.cdata[, -grep("_mut$", colnames(vali.cdata)), drop = FALSE],
  vali.cdata[, match(colnames(vali.mut.cdata), colnames(vali.cdata)), drop = FALSE],
  check.names = FALSE
)

vali.cdata$SurvObj <- Surv(
  vali.cdata$days,
  vali.cdata$os.status
)

ref.data <- vali.cdata
vali.mut.data <- vali.mut.cdata

if (nrow(ref.data) == 0) {
  stop("No samples remain after filtering.")
}

# ---- Cox model helper ----

get_cox_results <- function(age.levs, age.name) {
  
  data <- ref.data[
    ref.data$age_class_name %in% age.levs,
    ,
    drop = FALSE
  ]
  
  cat.var <- colnames(vali.mut.data)
  
  data$Sex <- data$Gender
  
  data <- data[
    !is.na(data$days) &
      data$First.Diagnosis %in% "Yes",
    ,
    drop = FALSE
  ]
  
  data <- na.omit(
    data[, c("SurvObj", cat.var, "age", "Sex", "Diag"), drop = FALSE]
  )
  
  if (nrow(data) == 0) {
    return(data.frame())
  }
  
  data$Male <- as.numeric(data$Sex %in% "Male")
  data$DMG <- as.numeric(data$Diag %in% "Midline")
  data$GBM <- as.numeric(data$Diag %in% "Glioblastoma")
  
  data.cat <- na.omit(
    data[, c(cat.var, "Male", "DMG", "GBM"), drop = FALSE]
  )
  
  data.cat.count <- apply(
    data.cat,
    2,
    function(x) sum(x %in% 1)
  )
  
  data.cat.count <- data.cat.count[data.cat.count >= 5]
  
  if (length(data.cat.count) == 0) {
    return(data.frame())
  }
  
  myformula <- as.formula(
    paste0(
      "SurvObj ~ ",
      paste(names(data.cat.count), collapse = " + "),
      " + age"
    )
  )
  
  myfit <- tryCatch(
    coxph(
      myformula,
      data = data,
      control = coxph.control(
        eps = 1e-09,
        toler.chol = .Machine$double.eps^0.75,
        iter.max = 20,
        toler.inf = sqrt(1e-09),
        outer.max = 50,
        timefix = TRUE
      )
    ),
    error = function(e) {
      message("Cox model failed for ", age.name, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(myfit)) {
    return(data.frame())
  }
  
  coef_tab <- data.frame(
    summary(myfit)$coef,
    check.names = FALSE
  )
  
  coef_tab$term <- row.names(coef_tab)
  row.names(coef_tab) <- NULL
  
  colnames(coef_tab) <- sub(
    "Pr\\(>\\|z\\|\\)",
    "p.value",
    colnames(coef_tab)
  )
  
  colnames(coef_tab) <- gsub(
    "exp\\(coef\\)",
    "exp.coef",
    colnames(coef_tab),
    fixed = TRUE
  )
  
  colnames(coef_tab) <- gsub(
    "se\\(coef\\)",
    "se.coef",
    colnames(coef_tab),
    fixed = TRUE
  )
  
  coef_tab$type <- gsub("_mut", "", coef_tab$term)
  coef_tab$age.name <- age.name
  coef_tab$n <- data.cat.count[match(coef_tab$term, names(data.cat.count))]
  coef_tab$n <- ifelse(is.na(coef_tab$n), nrow(data), coef_tab$n)
  
  coef_tab
}

# ---- Run age-stratified models ----

age_levels <- levels(ref.data$age_class_name)

cox_results_list <- list(
  get_cox_results(age.levs = age_levels[1], age.name = "PED"),
  get_cox_results(age.levs = age_levels[2], age.name = "ADO"),
  get_cox_results(age.levs = age_levels[1:2], age.name = "PED+ADO"),
  get_cox_results(age.levs = age_levels[3], age.name = "YA"),
  get_cox_results(age.levs = age_levels[2:3], age.name = "AYA"),
  get_cox_results(age.levs = age_levels[1:3], age.name = "PED+AYA"),
  get_cox_results(age.levs = age_levels[4], age.name = "ADULT"),
  get_cox_results(age.levs = age_levels[1:4], age.name = "ALL")
)

cox_results_list <- lapply(cox_results_list, function(x) {
  if (nrow(x) == 0) {
    return(x)
  }
  
  x$fdr.bonf <- p.adjust(
    x$p.value,
    method = "bonferroni"
  )
  
  x
})

cox_results <- do.call(rbind, cox_results_list)
cox_results <- data.frame(cox_results, check.names = FALSE)

if (nrow(cox_results) == 0) {
  stop("No Cox model results were produced.")
}

cox_results$signed.p <- -log10(cox_results$p.value) * sign(cox_results$coef)
cox_results$signed.fdr <- -log10(cox_results$fdr.bonf) * sign(cox_results$coef)

# ---- Format results for export and plotting ----

plotme <- data.frame(
  type = factor(
    cox_results$type,
    levels = c(
      gsub("_mut", "", colnames(vali.mut.data)),
      "age",
      "Male",
      "DMG",
      "GBM"
    )
  ),
  signed.p = cox_results$signed.p,
  signed.fdr = cox_results$signed.fdr,
  age.name = cox_results$age.name,
  n = cox_results$n,
  stringsAsFactors = FALSE
)

plotme <- plotme[
  !is.na(plotme$age.name) &
    !is.na(plotme$type),
  ,
  drop = FALSE
]

plotme <- plotme[
  order(abs(plotme$signed.p), decreasing = TRUE),
  ,
  drop = FALSE
]

plotme$type <- ifelse(
  plotme$type %in% "H33A",
  "H3-3A",
  as.character(plotme$type)
)

plotme$age.name <- factor(
  plotme$age.name,
  levels = c("PED", "ADO", "AYA", "PED+AYA", "ADULT", "ALL")
)

mylev <- unique(c(
  "GBM",
  "DMG",
  "age",
  "Male",
  as.character(unique(plotme[abs(plotme$signed.fdr) > 1, "type"]))
))

plotme$type <- factor(plotme$type, levels = rev(mylev))
plotme <- na.omit(plotme)

# ---- Save full results table ----

write.table(
  plotme,
  file = output_table,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# ---- Restrict plot to PED and AYA ----

plotme <- plotme[
  plotme$age.name %in% c("PED", "AYA"),
  ,
  drop = FALSE
]

if (nrow(plotme) == 0) {
  stop("No PED/AYA results available for plotting.")
}

# ---- Plot ----

p <- ggplot(
  plotme,
  aes(
    x = age.name,
    y = type,
    size = abs(signed.fdr),
    color = signed.fdr,
    label = n
  )
) +
  geom_point() +
  geom_point(
    data = subset(plotme, abs(signed.fdr) < 1),
    color = "gray80",
    fill = "white",
    pch = 21
  ) +
  geom_text(
    color = "black",
    size = 2,
    nudge_y = 0.5
  ) +
  scale_size(
    limits = c(0, max(abs(plotme$signed.fdr), na.rm = TRUE)),
    range = c(2, 8),
    breaks = c(0, 1, 2, 4, 6, 8, 10)
  ) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    legend.position = "right",
    legend.title.align = 0.5
  )

pdf(output_plot, height = 4, width = 3.5, useDingbats = FALSE)

print(
  p +
    scale_colour_gradientn(
      colours = c(
        muted("blue"),
        "white",
        "white",
        muted("red")
      ),
      values = scales::rescale(c(-5, -1, -0.2, 0.2, 1, 5)),
      limits = c(-5, 5),
      oob = squish,
      breaks = c(-5, -1, 0, 1, 5),
      labels = c("-5", "", "ns", "", "5"),
      name = "signed fdr",
      guide = guide_colorbar(
        direction = "vertical",
        ticks = TRUE,
        frame.colour = "black",
        frame.linewidth = 0.5,
        ticks.colour = "black",
        ticks.linewidth = 0.5
      )
    )
)

dev.off()

message("Wrote: ", output_table)
message("Wrote: ", output_plot)
