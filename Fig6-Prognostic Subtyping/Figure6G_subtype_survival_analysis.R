#!/usr/bin/env Rscript

# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

## =========================================================
## Figure 6G / STable 6 (Sheet 8: Subtype-Survival_Assoc)
##
## Reproduces subtype-associated survival analyses in the
## HOPE AYA HGG cohort.
##
## Data:
##   - ../data/STable1.xlsx, sheet ClinicalTable
##   - ../data/cDisc_mutation_10192023.tsv
##
## Outputs:
##   - Cox proportional hazards models (AG, CAG, CintAG)
##   - Kaplan–Meier curves with risk tables
##   - Exportable summary tables for manuscript
##
## Configurable:
##   - Subtype selection (cl.sel)
##   - Age cutoff (e.g., <=40, <=62)
##   - Inclusion of cancer group covariate
## =========================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(cowplot)
  library(colorspace)
  library(readxl)
})

clinical.data <- data.frame(
  read_xlsx(
    "../data/STable1.xlsx",
    sheet = "ClinicalTable",
    na = c("NA", "", "NaN")
  ),
  check.names = FALSE
)

mutation.raw <- read.delim(
  "../data/cDisc_mutation_10192023.tsv",
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "", "NaN")
)

mut.data <- mutation.raw[, -1, drop = FALSE]
row.names(mut.data) <- mutation.raw[, 1]


subtype.data <- data.frame(
  read_xlsx(
    "../data/STable6.xlsx",
    sheet = "Subtype-cDisc",
    na = c("NA", "", "NaN")
  ),
  check.names = FALSE
)


cluster.col <- c("1" = "#B94D3C", "2" = "#65B023", "3" = "#855F82")
male.cl.col <- darken(cluster.col, 0.2)
female.cl.col <- lighten(cluster.col, 0.2)
subtype.col <- c(male.cl.col, female.cl.col)

## --------------------------------------------
## prepare clinical data once
## --------------------------------------------
clinical.data <- as.data.frame(clinical.data)

clinical.data$cancer.group <- factor(clinical.data$Disc_Cancer_Group)
levels(factor(clinical.data$cancer.group))

map <- c("Other", "DMG", "HGG-NOS", "Other", "Other", "Other")
names(map) <- levels(factor(clinical.data$cancer.group))
levels(clinical.data$cancer.group) <- map
clinical.data$cancer.group <- relevel(clinical.data$cancer.group, ref = "Other")

## --------------------------------------------
## helper to extract Cox summary
## --------------------------------------------
extract_cox_summary <- function(model, model_name) {
  if (inherits(model, "coxph")) {
    summary_model <- summary(model)
    data.frame(
      Model = model_name,
      Variable = rownames(summary_model$coefficients),
      coef = summary_model$coefficients[, "coef"],
      HR = summary_model$conf.int[, "exp(coef)"],
      CI.low = summary_model$conf.int[, "lower .95"],
      CI.high = summary_model$conf.int[, "upper .95"],
      p.value = summary_model$coefficients[, "Pr(>|z|)"],
      check.names = FALSE,
      row.names = NULL
    )
  } else {
    data.frame(
      Model = model_name,
      Variable = "Model failed",
      coef = NA,
      HR = NA,
      CI.low = NA,
      CI.high = NA,
      p.value = NA,
      check.names = FALSE
    )
  }
}

## --------------------------------------------
## helper to format p-values for titles
## --------------------------------------------
format_p <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return("NA")
  }
  signif(x, digits)
}

## --------------------------------------------
## main function
## --------------------------------------------
run_survival_model <- function(
    clinical.data,
    mut.data,
    subtype.data,
    cl.sel = "2",
    age_cutoff = 40,
    include_cancer_group = FALSE,
    include_interaction_model = TRUE,
    max.day = 2000,
    text.size = 16
) {
  
  cl.sel <- as.character(cl.sel)
  data <- as.data.frame(clinical.data)
  
  ## ------------------
  ## derive variables
  ## ------------------
  data$H3F3A_mut <- as.matrix(mut.data)[
    row.names(mut.data) %in% "H3-3A",
    match(data$id, colnames(mut.data))
  ]
  
  data$is_initial <- as.numeric(
    data$cDisc_diagnosis_type %in% c("Initial CNS Tumor", "Progressive")
  )
  data$tumor.loc <- data$cDisc_tumor_loc
  data$Midline <- as.numeric(data$tumor.loc %in% "Midline")
  data$days <- suppressWarnings(as.numeric(data$cDisc_os))
  
  data$os_status_raw <- suppressWarnings(as.numeric(data$cDisc_os_status))
  data$os_days <- as.numeric(ifelse(data$days > max.day, max.day, data$days))
  data$os_status <- as.numeric(ifelse(data$days > max.day, 0, data$os_status_raw))
  data$os_status[is.na(data$days)] <- NA
  data$SurvObj <- with(data, Surv(os_days, os_status))
  
  data$subtype <- subtype.data[match(data$id, subtype.data$id), ]$protein.subtype
  data$cluster <- factor(gsub("(M|F)", "", data$subtype), levels = c("1", "2", "3"))
  data$age.class <- data$cDisc_age_class_name_derived
  data$ped <- as.numeric(data$age.class %in% "PED")
  data$Gender <- data$cDisc_Gender
  data$male <- as.numeric(data$Gender %in% "Male")
  data$cortical <- as.numeric(data$tumor.loc %in% "Cortical")
  data$dasc <- data$cDisc_clinical_status_at_collection_event
  data$age <- suppressWarnings(as.numeric(data$cDisc_age))
  data$treat_status <- data$cDisc_treat_status
  data$First.Diagnosis <- data$cDisc_First_Diagnosis
  data$diagnosis.type <- data$cDisc_diagnosis_type
  
  ## collapse selected cluster vs Other
  data$cl <- factor(
    ifelse(as.character(data$cluster) %in% cl.sel, cl.sel, "Other"),
    levels = c("Other", cl.sel)
  )
  
  ## ------------------
  ## subset
  ## ------------------
  data <- data[!is.na(data$days), ]
  data <- data[!data$dasc %in% "Deceased-due to disease", ]
  data <- data[data$treat_status %in% "Treatment naive", ]
  data <- data[data$First.Diagnosis %in% "Yes", ]
  data <- data[data$age <= age_cutoff, ]
  data <- data[data$age.class %in% c("PED", "ADO", "YA", "ADULT"), ]
  
  data$age.class <- factor(data$age.class)
  
  x <- table(data$cl)
  data <- data[data$cl %in% names(x[x >= 4]), ]
  data$cl <- droplevels(data$cl)
  
  if (!all(c("Other", cl.sel) %in% levels(data$cl))) {
    warning(paste0("After filtering, cl.sel = ", cl.sel, " does not have enough samples."))
    return(NULL)
  }
  
  data$age.class <- relevel(data$age.class, ref = "PED")
  data$primary <- as.numeric(data$diagnosis.type %in% "Primary")
  
  if (include_cancer_group) {
    data$cancer.group <- factor(data$cancer.group)
  }
  
  ## ------------------
  ## model formulas
  ## ------------------
  cox.ctrl <- coxph.control(
    eps = 1e-09,
    toler.chol = .Machine$double.eps^0.75,
    iter.max = 1000000,
    toler.inf = sqrt(1e-09),
    outer.max = 100,
    timefix = TRUE
  )
  
  if (include_cancer_group) {
    formulaAG <- as.formula(
      SurvObj ~ age.class + male + H3F3A_mut + Midline + is_initial + cancer.group
    )
    formulaCAG <- as.formula(
      SurvObj ~ cl + age.class + male + H3F3A_mut + Midline + is_initial + cancer.group
    )
  } else {
    formulaAG <- as.formula(
      SurvObj ~ age.class + male + H3F3A_mut + Midline + is_initial
    )
    formulaCAG <- as.formula(
      SurvObj ~ cl + age.class + male + H3F3A_mut + Midline + is_initial
    )
  }
  
  fitAG <- tryCatch(
    coxph(formulaAG, data = data, control = cox.ctrl),
    error = function(e) NULL
  )
  
  fitCAG <- tryCatch(
    coxph(formulaCAG, data = data, control = cox.ctrl),
    error = function(e) NULL
  )
  
  p.CAG_AG <- NA
  signed.p1 <- NA
  
  if (!is.null(fitCAG) && !is.null(fitAG)) {
    p.CAG_AG <- tryCatch(anova(fitCAG, fitAG)[2, 4], error = function(e) NA)
    if (!is.na(p.CAG_AG)) {
      signed.p1 <- log10(p.CAG_AG) * (-1) * sign(coef(summary(fitCAG))[1, 1])
    }
  }
  
  ## ------------------
  ## main adjusted curve
  ## ------------------
  new.df <- data.frame(
    cl = factor(levels(data$cl), levels = levels(data$cl)),
    age.class = factor(rep("PED", length(levels(data$cl))), levels = levels(data$age.class)),
    male = rep(0, length(levels(data$cl))),
    H3F3A_mut = rep(0, length(levels(data$cl))),
    Midline = rep(0, length(levels(data$cl))),
    is_initial = rep(0, length(levels(data$cl)))
  )
  
  if (include_cancer_group) {
    new.df$cancer.group <- factor(
      rep(levels(data$cancer.group)[1], nrow(new.df)),
      levels = levels(data$cancer.group)
    )
  }
  
  fit.curve <- tryCatch(
    survfit(fitCAG, newdata = new.df),
    error = function(e) NULL
  )
  
  splotCAG_AG <- NULL
  if (!is.null(fit.curve)) {
    splotCAG_AG <- suppressWarnings(
      suppressMessages(
        ggsurvplot(
          fit.curve,
          data = new.df,
          size = 0.5,
          conf.int = FALSE,
          risk.table = TRUE,
          legend = "none",
          title = paste0(
            "Model improvement P = ", format_p(p.CAG_AG),
            " | cl=", cl.sel,
            " | age<=", age_cutoff,
            ifelse(include_cancer_group, " | + cancer.group", "")
          ),
          ggtheme = theme_classic()
        )
      )
    )
  }
  
  ## ------------------
  ## optional interaction model
  ## ------------------
  fitCintAG <- NULL
  formulaCintAG <- NULL
  p.CintAG_AG <- NA
  p.interaction <- NA
  signed.p2 <- NA
  splotCintAG_AG <- NULL
  aligned <- NULL
  cox_results_CintG <- NULL
  
  if (include_interaction_model) {
    
    data$cl.male <- factor(
      ifelse(data$Gender %in% "Male", as.character(data$cl), "Other"),
      levels = c("Other", cl.sel)
    )
    
    if (include_cancer_group) {
      formulaCintAG <- as.formula(
        SurvObj ~ cl + cl.male + age.class + male + H3F3A_mut + Midline + is_initial + cancer.group
      )
    } else {
      formulaCintAG <- as.formula(
        SurvObj ~ cl + cl.male + age.class + male + H3F3A_mut + Midline + is_initial
      )
    }
    
    fitCintAG <- tryCatch(
      coxph(formulaCintAG, data = data, control = cox.ctrl),
      error = function(e) NULL
    )
    
    if (!is.null(fitCintAG) && !is.null(fitAG)) {
      p.CintAG_AG <- tryCatch(anova(fitCintAG, fitAG)[2, 4], error = function(e) NA)
      if (!is.na(p.CintAG_AG)) {
        signed.p2 <- log10(p.CintAG_AG) * (-1) * sign(sum(coef(summary(fitCintAG))[1:2, 1]))
      }
    }
    
    ## extract p-value for the sex-specific subtype effect term
    if (!is.null(fitCintAG)) {
      sm <- summary(fitCintAG)
      coef_names <- rownames(sm$coefficients)
      idx <- grep("^cl\\.male", coef_names)
      if (length(idx) > 0) {
        p.interaction <- sm$coefficients[idx[1], "Pr(>|z|)"]
      }
    }
    
    new.df2 <- data.frame(
      cl = factor(c(cl.sel, cl.sel, "Other"), levels = c("Other", cl.sel)),
      cl.male = factor(c(cl.sel, "Other", "Other"), levels = c("Other", cl.sel)),
      age.class = factor(rep("PED", 3), levels = levels(data$age.class)),
      male = c(0, 0, 0),
      H3F3A_mut = c(0, 0, 0),
      Midline = c(0, 0, 0),
      is_initial = c(0, 0, 0)
    )
    
    if (include_cancer_group) {
      new.df2$cancer.group <- factor(
        rep(levels(data$cancer.group)[1], nrow(new.df2)),
        levels = levels(data$cancer.group)
      )
    }
    
    fit.curve2 <- tryCatch(
      survfit(fitCintAG, newdata = new.df2),
      error = function(e) NULL
    )
    
    if (!is.null(fit.curve2)) {
      palette.use <- c(
        darken(cluster.col[[cl.sel]], 0.5),
        cluster.col[[cl.sel]],
        "gray"
      )
      
      splotCintAG_AG <- suppressWarnings(
        suppressMessages(
          ggsurvplot(
            fit.curve2,
            data = new.df2,
            size = 0.5,
            conf.int = FALSE,
            palette = palette.use,
            censor.size = 6,
            risk.table = TRUE,
            legend = "none",
            title = paste0(
              "Overall model improvement P = ", format_p(p.CintAG_AG),
              " | Sex-by-subtype interaction P = ", format_p(p.interaction),
              " | cl=", cl.sel,
              " | age<=", age_cutoff,
              ifelse(include_cancer_group, " | + cancer.group", "")
            ),
            ggtheme = theme_classic() + theme(
              plot.margin = margin(0, 0, 0, 0),
              plot.title = element_text(size = text.size),
              axis.title.x = element_text(size = text.size),
              axis.title.y = element_text(size = text.size),
              axis.text = element_text(size = text.size),
              legend.text = element_text(size = text.size),
              legend.title = element_text(size = text.size)
            ),
            risk.table.fontsize = 7,
            tables.theme = theme(
              plot.margin = margin(0, 0, 0, 0),
              text = element_text(size = text.size),
              axis.text = element_text(size = text.size),
              axis.title = element_text(size = text.size)
            )
          )
        )
      )
      
      aligned <- cowplot::align_plots(
        splotCintAG_AG$plot + theme(plot.margin = margin(r = 10)),
        splotCintAG_AG$table + theme(plot.margin = margin(r = 10)),
        align = "v"
      )
    }
    
    cox_results_CintG <- extract_cox_summary(
      fitCintAG,
      paste0(
        "Model CintAG | cl=", cl.sel,
        " | age<=", age_cutoff,
        ifelse(include_cancer_group, " | + cancer.group", "")
      )
    )
  }
  
  list(
    data = data,
    fitAG = fitAG,
    fitCAG = fitCAG,
    fitCintAG = fitCintAG,
    formulaAG = formulaAG,
    formulaCAG = formulaCAG,
    formulaCintAG = formulaCintAG,
    p.CAG_AG = p.CAG_AG,
    p.CintAG_AG = p.CintAG_AG,
    p.interaction = p.interaction,
    signed.p1 = signed.p1,
    signed.p2 = signed.p2,
    splotCAG_AG = splotCAG_AG,
    splotCintAG_AG = splotCintAG_AG,
    aligned = aligned,
    cox_results_CintG = cox_results_CintG
  )
}

## --------------------------------------------
## print helper
## --------------------------------------------
print_survival_result <- function(res, name = "model") {
  cat("\n============================\n")
  cat(name, "\n")
  cat("============================\n")
  
  if (is.null(res)) {
    cat("Result is NULL\n")
    return(invisible(NULL))
  }
  
  cat("formulaAG:\n")
  print(res$formulaAG)
  
  cat("\nformulaCAG:\n")
  print(res$formulaCAG)
  
  cat("\nfitCAG summary:\n")
  if (!is.null(res$fitCAG)) {
    print(summary(res$fitCAG))
  } else {
    cat("fitCAG is NULL\n")
  }
  
  cat("\nMain curve:\n")
  if (!is.null(res$splotCAG_AG)) {
    print(res$splotCAG_AG)
  } else {
    cat("splotCAG_AG is NULL\n")
  }
  
  cat("\nfitCintAG summary:\n")
  if (!is.null(res$fitCintAG)) {
    print(summary(res$fitCintAG))
  } else {
    cat("fitCintAG is NULL\n")
  }
  
  cat("\nInteraction p-value:\n")
  print(res$p.interaction)
  
  cat("\ncox_results_CintG:\n")
  if (!is.null(res$cox_results_CintG)) {
    print(res$cox_results_CintG)
  } else {
    cat("cox_results_CintG is NULL\n")
  }
  
  cat("\nInteraction curve:\n")
  if (!is.null(res$splotCintAG_AG)) {
    print(res$splotCintAG_AG)
  } else {
    cat("splotCintAG_AG is NULL\n")
  }
}

## --------------------------------------------
## save helper
## --------------------------------------------
save_survival_result <- function(res, prefix) {
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
script_dir <- resolve_script_dir("Figure6G_subtype_survival_analysis.R")
  output_dir <- file.path(script_dir, "output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(res$cox_results_CintG)) {
    write.table(
      res$cox_results_CintG,
      file.path(output_dir, paste0("TableS6_", prefix, "_cox_results_CintG.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  if (!is.null(res$splotCintAG_AG)) {
    pdf(file.path(output_dir, paste0("Figure6G_", prefix, "_interaction_survival.pdf")), width = 6, height = 7.5)
    print(res$splotCintAG_AG, newpage = FALSE)
    dev.off()
  }
}

## --------------------------------------------
## run models
## --------------------------------------------
res40_c1 <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "1",
  age_cutoff = 40,
  include_cancer_group = FALSE,
  include_interaction_model = TRUE
)

print_survival_result(res40_c1, "res40_c1")

res62_c2 <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "2",
  age_cutoff = 62,
  include_cancer_group = FALSE,
  include_interaction_model = TRUE
)

res40_c2 <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "2",
  age_cutoff = 40,
  include_cancer_group = FALSE,
  include_interaction_model = TRUE
)

res62_c2_cg <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "2",
  age_cutoff = 62,
  include_cancer_group = TRUE,
  include_interaction_model = TRUE
)

res40_c3 <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "3",
  age_cutoff = 40,
  include_cancer_group = FALSE,
  include_interaction_model = TRUE
)

res40_c2_cg <- run_survival_model(
  clinical.data = clinical.data,
  mut.data = mut.data,
  subtype.data = subtype.data,
  cl.sel = "2",
  age_cutoff = 40,
  include_cancer_group = TRUE,
  include_interaction_model = TRUE
)

## For supplementary Table 6: Subtype-Survival-Assoc
print_survival_result(res40_c1, "res40_c1")
print_survival_result(res62_c2, "res62_c2")
print_survival_result(res40_c2, "res40_c2")
print_survival_result(res40_c3, "res40_c3")
print_survival_result(res40_c2_cg, "res40_c2_cg")

## print Figure 6G
save_survival_result(res40_c2, "res40_c2")
