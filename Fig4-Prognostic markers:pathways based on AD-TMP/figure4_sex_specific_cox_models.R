#!/usr/bin/env Rscript

# ============================================================
# Figure 4 / Figure S4 | Sex-Specific Penalized Cox Models
# File: Figure4_sex_specific_cox_models.R
#
# Description:
#   Performs sex-stratified penalized Cox regression using glmnet
#   on proteogenomic data from pediatric and AYA high-grade glioma.
#
#   The script:
#     - Defines sex-specific candidate gene sets
#     - Fits ridge-penalized Cox models (alpha = 0)
#     - Uses fixed cross-validation folds for reproducibility
#     - Computes bootstrap confidence intervals for HRs and coefficients
#     - Derives combined prognostic scores
#     - Generates survival plots and hazard-ratio tables
#     - Creates figure-panel aliases after model execution
#
# Figures generated:
#   Main Figures:
#     - Figure 4E: Male model coefficients (weights plot)
#     - Figure 4I: PHKG2 survival plot (male)
#
#   Supplementary Figures:
#     - Figure S4D: Combined prognostic score survival (male)
#     - Figure S4E: Combined prognostic score survival (female)
#     - Figure S4F: Female model coefficients (weights plot)
#     - Figure S4I: MORC2 survival plot (male)
#
#   Supplementary Tables:
#     - STable4: Multivariable hazard ratio tables (male and female models)
#
# Inputs:
#   - ../data/STable1.xlsx, sheet ClinicalTable
#   - ../data/cDisc_mutation_10192023.tsv
#   - ../data/cDisc_proteome_imputed_data_09152023.tsv
#   - ../data/STable1.xlsx
#   - ../data/STable4.xlsx
#
# Outputs:
#   Written to the current figure folder.
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================



suppressPackageStartupMessages({
  library(glmnet)
  library(boot)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(readxl)
  library(scales)
})

if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA before running this script.", call. = FALSE)
}

set.seed(123)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

alpha_value <- 0
n_boot <- 10000
max_day <- 2000
data_type <- "cDiscovery"

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
script_dir <- resolve_script_dir("Figure4_sex_specific_cox_models.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- "../data"

sex_colors <- c(
  Male = "#0707CF",
  Female = "#CC0303",
  Both = "black"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

safe_scale_rows <- function(x) {
  x <- as.matrix(x)
  t(scale(t(x)))
}

save_pdf <- function(plot_object, filename, width, height) {
  pdf(file.path(output_dir, filename), width = width, height = height)
  print(plot_object)
  dev.off()
}

get_coef_plot_height <- function(n_predictors, min_height = 2.2, max_height = 8) {
  pmin(max_height, pmax(min_height, 0.35 * n_predictors + 1.2))
}

make_survival_data <- function(clinical, max_day = 2000) {
  os_time <- as.numeric(clinical$cDisc_os)
  pfs_time <- as.numeric(clinical$cDisc_pfs)
  os_event <- as.numeric(clinical$cDisc_os_status)
  pfs_event <- as.numeric(clinical$cDisc_pfs_status)

  os_status <- ifelse(os_time > max_day, 0, os_event)
  pfs_status <- ifelse(pfs_time > max_day, 0, pfs_event)
  
  os_days <- ifelse(os_time > max_day, max_day, os_time)
  pfs_days <- ifelse(pfs_time > max_day, max_day, pfs_time)
  
  data.frame(
    id = clinical$id,
    days = os_days,
    pfs.days = pfs_days,
    os.status = os_status,
    pfs.status = pfs_status,
    cohort = clinical$cohort,
    SurvObj = survival::Surv(os_days, os_status),
    PSurvObj = survival::Surv(pfs_days, pfs_status),
    stringsAsFactors = FALSE
  )
}

add_clinical_covariates <- function(surv_data, clinical) {
  clinical <- as.data.frame(clinical)
  rownames(clinical) <- clinical$id
  
  idx <- match(surv_data$id, rownames(clinical))
  
  surv_data %>%
    mutate(
      age = clinical$cDisc_age[idx],
      dasc = clinical$cDisc_clinical_status_at_collection_event[idx],
      treat_status = clinical$cDisc_treat_status[idx],
      First.Diagnosis = clinical$cDisc_First_Diagnosis[idx],
      age_class_new = clinical$cDisc_age_class_name_derived[idx],
      Gender = clinical$cDisc_Gender[idx],
      Cortical = clinical$cDisc_tumor_loc[idx] == "Cortical",
      Midline = as.numeric(clinical$cDisc_tumor_loc[idx] == "Midline"),
      Grade = as.numeric(clinical$cDisc_WHO_Grade[idx] %in% 3),
      tumor = as.numeric(
        clinical$cDisc_diagnosis_type[idx] %in% c("Primary", "Progressive")
      ),
      tumor = 1 - tumor,
      Male = as.numeric(Gender == "Male"),
      ado = as.numeric(age_class_new == "ADO")
    )
}

collapse_protein_to_gene <- function(protein_data, protein_anno) {
  gene_counts <- table(protein_anno$ApprovedGeneSymbol)
  
  protein_by_gene <- data.frame(
    do.call(
      rbind,
      lapply(names(gene_counts), function(gene) {
        apply(
          protein_data[protein_anno$ApprovedGeneSymbol == gene, , drop = FALSE],
          2,
          mean,
          na.rm = TRUE
        )
      })
    )
  )
  
  rownames(protein_by_gene) <- names(gene_counts)
  protein_by_gene
}

select_sex_signature <- function(stable4, sex_label) {
  sex_lower <- tolower(sex_label)
  
  cdisc_cols <- grep(
    paste0("comb.signed.log10p.cdisc.", sex_lower),
    colnames(stable4),
    value = TRUE
  )
  
  locfdr_cols <- grep(
    paste0("comb.signed.log10locfdr.ref.cont.", sex_lower),
    colnames(stable4),
    value = TRUE
  )
  
  sig_by <- apply(
    stable4[, cdisc_cols],
    2,
    function(x) {
      stable4$gene[p.adjust(10^((-1) * abs(x)), method = "BY") < 0.1]
    }
  )
  
  sig_by <- unique(as.character(unlist(sig_by[1:2])))
  
  cdisc_table <- stable4[, c("gene", cdisc_cols)]
  locfdr_table <- stable4[, c("gene", locfdr_cols)]
  
  selection_table <- merge(cdisc_table, locfdr_table, all = TRUE)
  
  ped_p <- paste0("PED.comb.signed.log10p.cdisc.", sex_lower)
  ped_l <- paste0("PED.comb.signed.log10locfdr.ref.cont.", sex_lower)
  ado_p <- paste0("ADO.comb.signed.log10p.cdisc.", sex_lower)
  ado_l <- paste0("ADO.comb.signed.log10locfdr.ref.cont.", sex_lower)
  
  sig_directional <- selection_table$gene[
    (
      sign(selection_table[[ped_l]]) == sign(selection_table[[ped_p]]) &
        abs(selection_table[[ped_p]]) > 2 &
        abs(selection_table[[ped_l]]) > 1
    ) |
      (
        sign(selection_table[[ado_l]]) == sign(selection_table[[ado_p]]) &
          abs(selection_table[[ado_p]]) > 2 &
          abs(selection_table[[ado_l]]) > 1
      )
  ]
  
  list(
    sig0 = sig_by,
    sig_directional = sig_directional,
    sig = unique(c(sig_directional, sig_by))
  )
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

clinical <- temporalCPSA::ageTMP_load_clinical(data_dir)
clinical$id <- temporalCPSA::ageTMP_normalize_sample_ids(clinical$id)
if (!"sample_id" %in% colnames(clinical)) {
  clinical$sample_id <- clinical$id
}
clinical$sample_id <- temporalCPSA::ageTMP_normalize_sample_ids(clinical$sample_id)

validation_clinical_raw <- read_xlsx(
  file.path(data_dir, "STable1.xlsx"),
  sheet = 5,
  na = "NA"
)

validation_protein_raw <- read_xlsx(
  file.path(data_dir, "STable1.xlsx"),
  sheet = 6,
  na = "NA"
)

stable4_sheet6 <- read_xlsx(
  file.path(data_dir, "STable4.xlsx"),
  sheet = 6,
  na = "NA"
)

# ------------------------------------------------------------------------------
# Define sex-specific candidate signatures
# ------------------------------------------------------------------------------

male_signature <- select_sex_signature(stable4_sheet6, "Male")
female_signature <- select_sex_signature(stable4_sheet6, "Female")

male_sig0 <- male_signature$sig0
male_sig_directional <- male_signature$sig_directional
male_sig <- male_signature$sig

female_sig0 <- female_signature$sig0
female_sig_directional <- female_signature$sig_directional
female_sig <- female_signature$sig

# Match overwritten legacy model definitions:
#   Male red   = c("IDH1", "IDH2", male_sig)
#   Female red = female_sig
male_candidate_genes <- unique(c("IDH1", "IDH2", na.omit(male_sig)))
female_candidate_genes <- unique(na.omit(female_sig))

both_candidate_genes <- unique(c(male_candidate_genes, female_candidate_genes))

both_mut_adj <- c("IDH1", "TP53", "ATRX", "H33A", "ATM")

# ------------------------------------------------------------------------------
# Prepare validation cohort data
# ------------------------------------------------------------------------------

validation_clinical <- data.frame(validation_clinical_raw)

validation_protein <- data.frame(validation_protein_raw)
rownames(validation_protein) <- validation_protein$Gene
validation_protein <- validation_protein[, colnames(validation_protein) != "Gene"]
colnames(validation_protein) <- gsub("\\.", "-", gsub("^X", "", colnames(validation_protein)))

validation_clinical <- validation_clinical %>%
  filter(
    age_class_derived %in% c("PED", "ADO"),
    Gender %in% c("Male", "Female")
  ) %>%
  mutate(
    tumor = 1 - is_initial_tumor,
    ado = as.numeric(age_class_derived == "ADO"),
    mut_H33A = H33A_mut
  ) %>%
  filter(surv_include)

validation_genes <- intersect(both_candidate_genes, rownames(validation_protein))

# ------------------------------------------------------------------------------
# Prepare cDiscovery survival data
# ------------------------------------------------------------------------------

rownames(clinical) <- clinical$id

data0 <- make_survival_data(clinical, max_day = max_day)
data0 <- add_clinical_covariates(data0, clinical)

survival_include_ids <- data0 %>%
  filter(
    !is.na(days),
    dasc != "Deceased-due to disease",
    treat_status == "Treatment naive",
    First.Diagnosis == "Yes"
  ) %>%
  pull(id)

data0 <- data0 %>%
  mutate(surv_include = id %in% survival_include_ids) %>%
  filter(
    surv_include,
    age_class_new %in% c("PED", "ADO"),
    Gender %in% c("Male", "Female")
  )

# ------------------------------------------------------------------------------
# Add mutation covariates
# ------------------------------------------------------------------------------

mutation_raw <- temporalCPSA::ageTMP_load_molecular(data_dir, modality = "mutation")
mutation_split <- temporalCPSA::ageTMP_split_annotation_matrix(
  mutation_raw,
  annotation_cols = 1,
  row_id = "ApprovedGeneSymbol"
)
mutation_data <- mutation_split$matrix
rownames(mutation_data) <- gsub("-", "", rownames(mutation_data))
mutation_data[is.na(mutation_data)] <- 0

available_mut_genes <- intersect(both_mut_adj, rownames(mutation_data))
missing_mut_genes <- setdiff(both_mut_adj, rownames(mutation_data))

if (length(missing_mut_genes) > 0) {
  message("Missing mutation genes: ", paste(missing_mut_genes, collapse = ", "))
}

sample_idx <- match(data0$id, colnames(mutation_data))
keep <- !is.na(sample_idx)

if (any(!keep)) {
  message("Dropping ", sum(!keep), " samples not found in mutation matrix")
}

mutation_sub <- mutation_data[
  available_mut_genes,
  sample_idx[keep],
  drop = FALSE
]

mutation_sub <- t(mutation_sub)

rownames(mutation_sub) <- data0$id[keep]
colnames(mutation_sub) <- paste0("mut_", available_mut_genes)

data0 <- data0[keep, , drop = FALSE]
mutation_sub[is.na(mutation_sub)] <- 0

data0 <- cbind(data0, mutation_sub)

stopifnot(nrow(data0) == nrow(mutation_sub))
stopifnot(all(rownames(mutation_sub) == data0$id))

# ------------------------------------------------------------------------------
# Add protein covariates
# ------------------------------------------------------------------------------

protein_raw <- temporalCPSA::ageTMP_load_molecular(data_dir, modality = "protein")
protein_split <- temporalCPSA::ageTMP_split_annotation_matrix(
  protein_raw,
  annotation_cols = 1:4,
  row_id = "ApprovedGeneSymbol"
)
protein_data <- protein_split$matrix
protein_anno <- protein_split$annotation

protein_by_gene <- collapse_protein_to_gene(protein_data, protein_anno)

protein_by_gene <- protein_by_gene[
  rowSums(is.na(protein_by_gene)) == 0,
  ,
  drop = FALSE
]

protein_sample_idx <- match(
  temporalCPSA::ageTMP_normalize_sample_ids(colnames(protein_by_gene)),
  clinical$id
)
if (any(is.na(protein_sample_idx))) {
  protein_sample_idx[is.na(protein_sample_idx)] <- match(
    temporalCPSA::ageTMP_normalize_sample_ids(colnames(protein_by_gene)[is.na(protein_sample_idx)]),
    clinical$sample_id
  )
}
colnames(protein_by_gene) <- clinical$id[protein_sample_idx]
protein_by_gene <- protein_by_gene[, !is.na(colnames(protein_by_gene)), drop = FALSE]

protein_genes <- intersect(both_candidate_genes, rownames(protein_by_gene))

protein_sub <- protein_by_gene[
  protein_genes,
  data0$id,
  drop = FALSE
]

protein_sub <- safe_scale_rows(protein_sub)
protein_sub <- t(protein_sub)

protein_interactions <- apply(
  protein_sub,
  2,
  function(x) ifelse(data0$ado == 1, x, 0)
)

colnames(protein_interactions) <- paste0("ado.", colnames(protein_interactions))

data0 <- data.frame(data0, protein_sub, protein_interactions)

male_data0 <- data0[data0$Gender == "Male", ]
female_data0 <- data0[data0$Gender == "Female", ]

# ------------------------------------------------------------------------------
# Sex-specific Cox model function
# ------------------------------------------------------------------------------

run_sex_specific_model <- function(sex_selected) {
  
  message("Running model for: ", sex_selected)
  
  if (sex_selected == "Male") {
    model_data <- male_data0
    candidate_genes <- male_candidate_genes
    sex_sig <- male_sig
    sex_sig_directional <- male_sig_directional
  } else if (sex_selected == "Female") {
    model_data <- female_data0
    candidate_genes <- female_candidate_genes
    sex_sig <- female_sig
    sex_sig_directional <- female_sig_directional
  } else {
    stop("sex_selected must be 'Male' or 'Female'.")
  }
  
  model_genes <- intersect(candidate_genes, validation_genes)
  model_genes <- intersect(model_genes, colnames(model_data))
  
  mutation_covariates <- "mut_H33A"
  clinical_covariates <- c("ado", "Midline", "Grade", "tumor")
  
  model_covariates <- c(
    mutation_covariates,
    model_genes,
    clinical_covariates
  )
  
  model_covariates <- intersect(model_covariates, colnames(model_data))
  
  # Only require complete clinical/mutation covariates.
  # Do NOT filter on all protein predictors, or the cohort collapses.
  base_covariates <- intersect(
    c("mut_H33A", "ado", "Midline", "Grade", "tumor"),
    colnames(model_data)
  )
  
  message(
    sex_selected,
    " N before filtering: ",
    nrow(model_data)
  )
  
  model_data <- model_data[
    complete.cases(model_data[, base_covariates, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  message(
    sex_selected,
    " N after clinical/mutation filtering: ",
    nrow(model_data)
  )
  
  X <- scale(as.matrix(model_data[, model_covariates, drop = FALSE]))
  Y <- model_data$SurvObj
  
  set.seed(123)
  foldid <- sample(rep(seq_len(10), length.out = nrow(X)))
  
  cv_fit <- cv.glmnet(
    X,
    Y,
    family = "cox",
    alpha = alpha_value,
    foldid = foldid
  )
  
  best_lambda <- cv_fit$lambda.min
  lambda_1se <- cv_fit$lambda.1se
  
  message(sex_selected, " optimal lambda: ", best_lambda)
  message(sex_selected, " lambda within 1 SE: ", lambda_1se)
  
  final_fit <- glmnet(
    X,
    Y,
    family = "cox",
    alpha = alpha_value,
    lambda = best_lambda
  )
  
  model_coefficients <- as.matrix(coef(final_fit))
  
  selected_variables <- rownames(model_coefficients)[
    abs(model_coefficients[, 1]) > 0
  ]
  
  # --------------------------------------------------------------------------
  # Bootstrap HRs
  # --------------------------------------------------------------------------
  
  bootstrap_hr <- function(data, indices) {
    X_boot <- X[indices, , drop = FALSE]
    Y_boot <- Y[indices]
    
    model <- glmnet(
      X_boot,
      Y_boot,
      family = "cox",
      alpha = alpha_value,
      lambda = best_lambda
    )
    
    exp(as.matrix(coef(model))[, 1])
  }
  
  set.seed(123)
  
  hr_boot <- boot(
    data = seq_len(nrow(X)),
    statistic = bootstrap_hr,
    R = n_boot
  )
  
  hr_data <- data.frame(
    Predictor = colnames(X),
    HR = apply(hr_boot$t, 2, median, na.rm = TRUE),
    CI_Lower = apply(hr_boot$t, 2, quantile, probs = 0.05, na.rm = TRUE),
    CI_Upper = apply(hr_boot$t, 2, quantile, probs = 0.95, na.rm = TRUE),
    n = apply(!is.na(hr_boot$t), 2, sum),
    sel = apply(hr_boot$t, 2, function(x) mean(x != 0, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  rownames(hr_data) <- hr_data$Predictor
  
  if (alpha_value == 0) {
    selected_from_ci <- hr_data$Predictor[
      hr_data$CI_Lower > 1 | hr_data$CI_Upper < 1
    ]
    
    selected_genes <- intersect(selected_from_ci, model_genes)
  } else {
    selected_genes <- intersect(selected_variables, model_genes)
  }
  
  message(sex_selected, " selected genes: ", paste(selected_genes, collapse = ", "))
  
  hr_data <- hr_data %>%
    arrange(HR) %>%
    mutate(Predictor = factor(Predictor, levels = unique(Predictor)))
  
  hr_plot <- ggplot(
    hr_data[hr_data$Predictor %in% selected_genes, ],
    aes(x = Predictor, y = HR)
  ) +
    geom_point(color = sex_colors[[sex_selected]], size = 3) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper),
      color = sex_colors[[sex_selected]],
      width = 0.2
    ) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    coord_flip() +
    scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
    labs(
      title = paste0("Hazard Ratios with 90% Bootstrap CIs (", sex_selected, ")"),
      x = "Predictors",
      y = "HR"
    ) +
    theme_bw() +
    theme(plot.title = element_text(size = 8))
 
  # --------------------------------------------------------------------------
  # Bootstrap coefficients
  # --------------------------------------------------------------------------
  
  bootstrap_coef <- function(data, indices) {
    X_boot <- X[indices, , drop = FALSE]
    Y_boot <- Y[indices]
    
    model <- glmnet(
      X_boot,
      Y_boot,
      family = "cox",
      alpha = alpha_value,
      lambda = best_lambda
    )
    
    as.matrix(coef(model))[, 1]
  }
  
  set.seed(123)
  
  coef_boot <- boot(
    data = seq_len(nrow(X)),
    statistic = bootstrap_coef,
    R = n_boot
  )
  
  coef_data <- data.frame(
    Predictor = colnames(X),
    coef = as.matrix(coef(final_fit))[, 1],
    CI_Lower = apply(coef_boot$t, 2, quantile, probs = 0.05, na.rm = TRUE),
    CI_Upper = apply(coef_boot$t, 2, quantile, probs = 0.95, na.rm = TRUE),
    n = apply(!is.na(coef_boot$t), 2, sum),
    sel = apply(coef_boot$t, 2, function(x) mean(x != 0, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  rownames(coef_data) <- coef_data$Predictor
  
  coef_data <- coef_data %>%
    filter(Predictor %in% selected_genes) %>%
    arrange(coef) %>%
    mutate(Predictor = factor(Predictor, levels = unique(Predictor)))
  
  write.table(
    coef_data,
    file = file.path(output_dir, paste0("coef_data_", sex_selected, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  coef_plot <- ggplot(coef_data, aes(y = Predictor, x = coef)) +
    geom_col(fill = sex_colors[[sex_selected]], width = 0.7) +
    geom_errorbar(
      aes(xmin = CI_Lower, xmax = CI_Upper),
      color = sex_colors[[sex_selected]],
      width = 0.2
    ) +
    geom_vline(xintercept = 0, linetype = 2, color = "black") +
    coord_cartesian(xlim = c(-0.15, 0.15)) +
    labs(
      title = paste0("Coefficients with 90% Bootstrap CIs (", sex_selected, ")"),
      x = "Coefficient",
      y = "Predictor"
    ) +
    theme_gray() +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  save_pdf(
    coef_plot,
    filename = paste0("weights_sel_", data_type, "_", sex_selected, ".pdf"),
    height = get_coef_plot_height(nrow(coef_data)),
    width = 6.5
  )
  
  selected_coefficients <- model_coefficients[
    rownames(model_coefficients) %in% selected_genes,
    ,
    drop = FALSE
  ]
  
  if (sex_selected == "Male") {
    male_coefficients <- selected_coefficients
    save(
      male_sig,
      male_sig_directional,
      male_coefficients,
      file = file.path(output_dir, "male_coefficients.RData")
    )
  }
  
  if (sex_selected == "Female") {
    female_coefficients <- selected_coefficients
    save(
      female_sig,
      female_sig_directional,
      female_coefficients,
      file = file.path(output_dir, "female_coefficients.RData")
    )
  }
  
  # --------------------------------------------------------------------------
  # Combined prognostic score
  # --------------------------------------------------------------------------
  
  score_genes <- intersect(selected_genes, model_genes)
  
  if (length(score_genes) > 0) {
    
    score_weights <- as.numeric(model_coefficients[score_genes, 1])
    names(score_weights) <- score_genes
    
    model_data$prog.combined <- as.numeric(
      scale(as.matrix(model_data[, score_genes, drop = FALSE]) %*% score_weights * -1)
    )
    
    model_data$cl <- cut(
      model_data$prog.combined,
      breaks = c(
        min(model_data$prog.combined, na.rm = TRUE),
        median(model_data$prog.combined, na.rm = TRUE),
        max(model_data$prog.combined, na.rm = TRUE)
      ),
      include.lowest = TRUE,
      labels = c("low", "high")
    )
    
    adjustment_covariates <- model_covariates[
      !model_covariates %in% model_genes
    ]
    
    score_fit <- coxph(
      as.formula(
        paste0(
          "SurvObj ~ prog.combined + ",
          paste(adjustment_covariates, collapse = " + ")
        )
      ),
      data = model_data
    )
    
    score_ci <- exp(confint(score_fit, level = 0.90))
    
    plot_data <- data.frame(
      Variable = rownames(score_ci),
      HR = exp(coef(score_fit)),
      CI_lower = score_ci[, 1],
      CI_upper = score_ci[, 2],
      stringsAsFactors = FALSE
    )
    
    plot_data <- plot_data %>%
      arrange(HR) %>%
      mutate(Variable = factor(Variable, levels = unique(Variable)))
    
    write.table(
      plot_data,
      file = file.path(
        output_dir,
        paste0("multivariate_hr_table_", sex_selected, "_", data_type, ".tsv")
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    hr_plot <- ggplot(plot_data, aes(x = Variable, y = HR)) +
      geom_point(color = sex_colors[[sex_selected]], size = 3) +
      geom_errorbar(
        aes(ymin = CI_lower, ymax = CI_upper),
        width = 0.2,
        color = sex_colors[[sex_selected]]
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      scale_y_log10() +
      labs(
        x = "Variables",
        y = "Hazard Ratio"
      ) +
      facet_wrap(~paste0("HRs with 90% CIs (", data_type, ")")) +
      theme_bw()
    
    save_pdf(
      hr_plot,
      filename = paste0("multivariate_hr_plot_", sex_selected, "_", data_type, ".pdf"),
      width = 4,
      height = 3
    )
    
    cluster_fit <- coxph(
      as.formula(
        paste0(
          "SurvObj ~ cl + ",
          paste(adjustment_covariates, collapse = " + ")
        )
      ),
      control = coxph.control(iter.max = 2000),
      data = model_data
    )
    
    new_data <- data.frame(
      cl = factor(c("high", "low"), levels = c("low", "high"))
    )
    
    rownames(new_data) <- paste0(c("high", "low"), " N = ", table(model_data$cl))
    
    adjustment_data <- data.frame(
      matrix(
        0,
        nrow = 2,
        ncol = length(adjustment_covariates),
        dimnames = list(NULL, adjustment_covariates)
      )
    )
    
    new_data <- data.frame(new_data, adjustment_data)
    
    combined_survfit <- survfit(cluster_fit, newdata = new_data)
    
    combined_surv_plot <- ggsurvplot(
      combined_survfit,
      data = new_data,
      size = 0.5,
      title = paste0("Combined score - ", data_type, " (", sex_selected, ")"),
      conf.int = FALSE,
      linetype = "strata",
      pval = paste0(
        "pval (cont): ",
        round(summary(score_fit)$coef["prog.combined", 5], 3)
      ),
      palette = c(sex_colors[[sex_selected]], "gray40"),
      risk.table = TRUE,
      risk.table.fontsize = 5,
      ggtheme = theme_bw(),
      tables.theme = theme_bw()
    )
    
    aligned_combined <- cowplot::align_plots(
      combined_surv_plot$plot + theme(plot.margin = margin(r = 10)),
      combined_surv_plot$table + theme(plot.margin = margin(r = 10)),
      align = "v"
    )
    
    pdf(
      file.path(output_dir, paste0("prog_combined_", sex_selected, "_", data_type, ".pdf")),
      width = 6.2,
      height = 6.8
    )
    
    print(
      cowplot::plot_grid(
        aligned_combined[[1]],
        aligned_combined[[2]],
        ncol = 1,
        rel_heights = c(4, 1)
      )
    )
    
    dev.off()
  }
  
  # --------------------------------------------------------------------------
  # Male-only PHKG2/MORC2 survival plots
  # --------------------------------------------------------------------------
  
  if (sex_selected == "Male") {
    
    predictors_to_plot <- intersect(c("PHKG2", "MORC2"), colnames(model_data))
    
    plot_single_predictor <- function(pred) {
      
      adjustment_covs <- model_covariates[
        !model_covariates %in% pred &
          !model_covariates %in% model_genes
      ]
      
      continuous_fit <- coxph(
        as.formula(
          paste0("SurvObj ~ ", pred, " + ", paste(adjustment_covs, collapse = " + "))
        ),
        data = model_data
      )
      
      null_fit <- coxph(
        as.formula(
          paste0("SurvObj ~ ", paste(adjustment_covs, collapse = " + "))
        ),
        data = model_data
      )
      
      temp_data <- model_data
      
      temp_data$cl <- factor(
        cut(
          temp_data[, pred],
          breaks = c(
            min(temp_data[, pred], na.rm = TRUE),
            median(temp_data[, pred], na.rm = TRUE),
            max(temp_data[, pred], na.rm = TRUE)
          ),
          include.lowest = TRUE,
          labels = c("low", "high")
        ),
        levels = c("low", "high")
      )
      
      anova_p <- anova(continuous_fit, null_fit)[2, 4]
      
      binary_fit <- coxph(
        as.formula(
          paste0("SurvObj ~ cl + ", paste(adjustment_covs, collapse = " + "))
        ),
        data = temp_data
      )
      
      new_data <- data.frame(
        cl = factor(c("high", "low"), levels = c("low", "high"))
      )
      
      rownames(new_data) <- paste0(c("high", "low"), " N = ", table(temp_data$cl))
      
      adjustment_data <- data.frame(
        matrix(
          0,
          nrow = 2,
          ncol = length(adjustment_covs),
          dimnames = list(NULL, adjustment_covs)
        )
      )
      
      new_data <- data.frame(new_data, adjustment_data)
      
      predictor_survfit <- survfit(binary_fit, newdata = new_data)
      
      predictor_plot <- ggsurvplot(
        predictor_survfit,
        data = new_data,
        title = paste0(pred, " - ", data_type, " (", sex_selected, ")"),
        conf.int = FALSE,
        linetype = "strata",
        pval = paste0(
          "pval (cont) = ",
          round(summary(continuous_fit)$coef[pred, 5], 3)
        ),
        palette = c(sex_colors[[sex_selected]], "gray40"),
        risk.table = TRUE,
        risk.table.fontsize = 5,
        size = .5,
        ggtheme = theme_bw(),
        tables.theme = theme_bw()
      )
      aligned_predictor <- cowplot::align_plots(
        predictor_plot$plot + theme(plot.margin = margin(r = 10)),
        predictor_plot$table + theme(plot.margin = margin(r = 10)),
        align = "v"
      )
      
      output_file <- paste0(pred, "_", sex_selected, "_", data_type, "_splot.pdf")
      
      pdf(file.path(output_dir, output_file), width = 6.2, height = 6.8)
      
      print(
        cowplot::plot_grid(
          aligned_predictor[[1]],
          aligned_predictor[[2]],
          ncol = 1,
          rel_heights = c(4, 1)
        )
      )
      
      dev.off()
    }
    
    lapply(predictors_to_plot, plot_single_predictor)
  }
  
  invisible(list(
    sex = sex_selected,
    model_ids = model_data$id,
    model_covariates = model_covariates,
    foldid = foldid,
    best_lambda = best_lambda,
    selected_genes = selected_genes,
    coef_data = coef_data,
    hr_data = hr_data,
    coefficients = selected_coefficients
  ))
}

# ------------------------------------------------------------------------------
# Run models
# ------------------------------------------------------------------------------

sex_model_results <- lapply(c("Male", "Female"), run_sex_specific_model)
names(sex_model_results) <- c("Male", "Female")

# ------------------------------------------------------------------------------
# Figure panel aliases
# ------------------------------------------------------------------------------

figure_aliases <- c(
  "weights_sel_cDiscovery_Male.pdf" =
    "Figure4E_weights_sel_cDiscovery_Male.pdf",
  "prog_combined_Male_cDiscovery.pdf" =
    "FigureS4D_prog_combined_Male_cDiscovery.pdf",
  "PHKG2_Male_cDiscovery_splot.pdf" =
    "Figure4I_PHKG2_Male_cDiscovery_splot.pdf",
  "MORC2_Male_cDiscovery_splot.pdf" =
    "FigureS4I_MORC2_Male_cDiscovery_splot.pdf",

  "multivariate_hr_table_Male_cDiscovery.tsv" =
    "STable4_multivariate_hr_table_Male_cDiscovery.tsv",
  
  "weights_sel_cDiscovery_Female.pdf" =
    "FigureS4F_weights_sel_cDiscovery_Female.pdf",
  "prog_combined_Female_cDiscovery.pdf" =
    "FigureS4E_prog_combined_Female_cDiscovery.pdf",

  "multivariate_hr_table_Female_cDiscovery.tsv" =
    "STable4_multivariate_hr_table_Female_cDiscovery.tsv"
)

for (source_file in names(figure_aliases)) {
  source_path <- file.path(output_dir, source_file)
  target_path <- file.path(output_dir, figure_aliases[[source_file]])
  
  if (file.exists(source_path)) {
    file.copy(source_path, target_path, overwrite = TRUE)
  } else {
    warning("Expected source file not found: ", source_path)
  }
}

message("Done. Male and Female outputs written to: ", output_dir)


analysis_ready_data <- list(
  
  cDiscovery = list(
    data = data0,
    male_data = male_data0,
    female_data = female_data0
  ),
  
  validation = list(
    clinical = validation_clinical,
    protein = validation_protein,
    
    # Explicit downstream-ready survival object
    validation_data = validation_clinical %>%
      mutate(
        SurvObj = survival::Surv(survival_days, os_status),
        tumor = 1 - is_initial_tumor,
        ado = as.numeric(age_class_derived == "ADO"),
        Male = as.numeric(Gender == "Male"),
        mut_H33A = H33A_mut
      )
  ),
  
  signatures = list(
    male_sig = male_sig,
    female_sig = female_sig,
    male_candidate_genes = male_candidate_genes,
    female_candidate_genes = female_candidate_genes,
    both_candidate_genes = both_candidate_genes,
    validation_genes = validation_genes
  ),
  
  coefficients = list(
    male_coefficients = sex_model_results$Male$coefficients,
    female_coefficients = sex_model_results$Female$coefficients
  ),
  
  parameters = list(
    max_day = max_day,
    alpha_value = alpha_value,
    data_type = data_type
  )
)

saveRDS(
  analysis_ready_data,
  file = file.path(output_dir, "Figure4_analysis_ready_data.rds")
)

message(
  "Saved downstream analysis-ready data to: ",
  file.path(output_dir, "Figure4_analysis_ready_data.rds")
)
