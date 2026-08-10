#!/usr/bin/env Rscript

# ============================================================
# Figure 1G | Reference Cohort Survival Analysis by Age Class
# File: Figure1G_reference_survival_assoc.R
#
# Description:
#   Generates Kaplan–Meier survival curves for the external
#   reference cohort stratified by age class and exports a
#   supplementary survival summary table.
#
# Inputs:
#   - ../data/STable1.xlsx, sheets ClinicalTable and Ref_ClinicalTable
#
# Outputs:
#   - Figure1G_reference_survival_assoc.pdf
#   - Both_supplementary_survival_table.csv
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(readxl)
  library(RColorBrewer)
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
script_dir <- resolve_script_dir("Figure1G_reference_survival_assoc.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

get_green <- function(n) {
  colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-(1:2)])(n)
}


clinical.data <- data.frame(
  read_xlsx(
    "../data/STable1.xlsx",
    sheet = "ClinicalTable",
    na = c("NA", "", "NaN")
  ),
  check.names = FALSE
)
row.names(clinical.data)=clinical.data$id



ref.data <- data.frame(
  read_xlsx(
    "../data/STable1.xlsx",
    sheet = "Ref_ClinicalTable",
    na = c("NA", "", "NaN")
  ),
  check.names = FALSE
)
ref.data$age_class_name=factor(ref.data$age_class_name,levels=c("PED","ADO","YA","ADULT","SEN"))

ref.data$SurvObj=Surv(
  ref.data$days,
  ref.data$os.status
)
vali.surv.data <- ref.data[, c("id", "SurvObj", "age_class_name", "Gender")]
vali.surv.data <- data.frame(vali.surv.data, as.matrix(vali.surv.data$SurvObj))


surv.data <- vali.surv.data
surv.data$SurvObj <- with(surv.data, Surv(time, status))

get_surv.cl <- function(
    cl = NULL,
    mytitle = "Age class",
    data = surv.data,
    scheme = get_green,
    mypal0 = NULL,
    cl.name = "age_class_name"
) {
  data <- data[!is.na(data$time), ]
  
  data$cl <- data[, cl.name]
  
  fit <- survfit(SurvObj ~ cl, data = data)
  fit.table <- as.data.frame(summary(fit)$table, check.names = FALSE)
  
  surv_table <- data.frame(
    Group = row.names(fit.table),
    n = fit.table[, "records"],
    RMST = round(fit.table[, "rmean"], 1),
    se = round(fit.table[, "se(rmean)"], 2),
    Events = fit.table[, "events"],
    Censored = fit.table[, "n.start"] - fit.table[, "events"]
  )
  
  log_rank_test <- survdiff(Surv(time, status) ~ cl, data = data)
  log_p_value <- pchisq(
    log_rank_test$chisq,
    length(log_rank_test$n) - 1,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  surv_table <- surv_table %>%
    dplyr::mutate(
      log_p_value = log_p_value,
      RMST_SE = paste0(RMST, " (", se, ")")
    ) %>%
    dplyr::select(Group, n, RMST_SE, Events, Censored, log_p_value)
  
  print(surv_table)
  
  write.csv(
    surv_table,
    file.path(output_dir, paste0(mytitle, "_supplementary_survival_table.csv")),
    row.names = FALSE
  )
  
  if (!is.null(scheme)) {
    mypal0 <- scheme(length(levels(factor(data$cl))))
    names(mypal0) <- levels(data$cl)
  }
  
  mypal <- mypal0[match(levels(data$cl), names(mypal0))]
  
  splot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    palette = as.character(mypal),
    conf.int = FALSE,
    risk.table = TRUE,
    size = 0.5,
    title = mytitle
  )
  
  splot
}

surv.data$Gender <- factor(surv.data$Gender, levels = c("Male", "Female"))
surv.data <- surv.data[!surv.data$id %in% clinical.data$id, ]

splot <- get_surv.cl(mytitle = "Both")

aligned <- cowplot::align_plots(
  splot$plot + theme(plot.margin = margin(r = 10)),
  splot$table + theme(plot.margin = margin(r = 10)),
  align = "v"
)

pdf(file.path(output_dir, "Figure1G_reference_survival_assoc.pdf"), width = 6, height = 7)
print(cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = c(4, 1.5)))
dev.off()
