#!/usr/bin/env Rscript

# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

# Figure 1F manuscript reproduction. The age-class structure shown in this
# panel was derived early in the study from the data versions available during
# manuscript development. This script preserves the archived AD-TMP heatmap
# input and the descriptive survival-days annotation used to render the
# published visualization. For a companion analysis based on the final
# repository data tables, run Figure1F_AD_TMP_heatmap_source_derived.R.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]]))
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop("Cannot determine script directory for `", script_name, "`. Run from the script folder or with Rscript.", call. = FALSE)
}
script_dir <- resolve_script_dir("Figure1F_AD_TMP_heatmap.R")

data_dir <- Sys.getenv("AGETMP_DATA_DIR", file.path(repo_root, "data"))
out_dir <- Sys.getenv("AGETMP_OUTPUT_DIR", file.path(script_dir, "output"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

figure1f_mode <- Sys.getenv("AGETMP_FIGURE1F_MODE", "manuscript")
if (!figure1f_mode %in% c("manuscript", "source")) {
  stop("`AGETMP_FIGURE1F_MODE` must be `manuscript` or `source`.", call. = FALSE)
}
output_prefix <- if (identical(figure1f_mode, "manuscript")) {
  "figure1f_ad_tmp_heatmap_manuscript"
} else {
  "figure1f_ad_tmp_heatmap_source_derived"
}

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists(file.path(repo_root, "temporalCPSA"))) {
  pkgload::load_all(file.path(repo_root, "temporalCPSA"), quiet = TRUE)
} else if (requireNamespace("temporalCPSA", quietly = TRUE)) {
  library(temporalCPSA)
} else {
  stop("Install temporalCPSA or run this script from a repository root containing temporalCPSA/.", call. = FALSE)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) default else out
}

fill_edge_na_by_row <- function(mat) {
  mat <- as.matrix(mat)
  for (i in seq_len(nrow(mat))) {
    x <- mat[i, ]
    ok <- which(!is.na(x))
    if (length(ok) == 0) {
      next
    }
    x[is.na(x)] <- x[ok[[1]]]
    mat[i, ] <- x
  }
  mat
}

read_feature_matrix_file <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE)
  rownames(x) <- x[[1]]
  x <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(x) <- "numeric"
  x
}

generate_ad_tmp_grid <- function(data_dir, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fit_age_max <- env_num("AGETMP_FIT_AGE_MAX", 80)
  center_age_max <- env_num("AGETMP_CENTER_AGE_MAX", 80)
  span <- env_num("AGETMP_TRAJECTORY_SPAN", 1.5)
  span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_TRAJECTORY_SPAN", span),
    rna = env_num("AGETMP_RNA_TRAJECTORY_SPAN", span),
    phospho = env_num("AGETMP_PHOSPHO_TRAJECTORY_SPAN", span)
  )
  adaptive_span <- !identical(Sys.getenv("AGETMP_ADAPTIVE_TRAJECTORY_SPAN", unset = "FALSE"), "FALSE")
  min_span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_MIN_TRAJECTORY_SPAN", 0.5),
    rna = env_num("AGETMP_RNA_MIN_TRAJECTORY_SPAN", 1.0),
    phospho = env_num("AGETMP_PHOSPHO_MIN_TRAJECTORY_SPAN", 0.5)
  )
  max_span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_MAX_TRAJECTORY_SPAN", 3.0),
    rna = env_num("AGETMP_RNA_MAX_TRAJECTORY_SPAN", 3.0),
    phospho = env_num("AGETMP_PHOSPHO_MAX_TRAJECTORY_SPAN", 3.0)
  )
  span_step <- env_num("AGETMP_TRAJECTORY_SPAN_STEP", 0.1)
  feature_limit <- env_num("AGETMP_FEATURE_LIMIT", Inf)
  dynamic_proportion <- env_num("AGETMP_DYNAMIC_PROPORTION", 0.5)
  max_dynamic_features <- env_num("AGETMP_MAX_DYNAMIC_FEATURES", 5000)
  figure1f_max_features_by_modality <- c(
    protein = env_num("AGETMP_FIGURE1F_PROTEIN_FEATURES", 5000),
    phospho = env_num("AGETMP_FIGURE1F_PHOSPHO_FEATURES", 3804),
    rna = env_num("AGETMP_FIGURE1F_RNA_FEATURES", 5000)
  )
  grid_step <- env_num("AGETMP_GRID_STEP", 0.5)
  n_cores <- as.integer(env_num("AGETMP_N_CORES", 1))
  remove_dasc <- !identical(Sys.getenv("AGETMP_REMOVE_DASC", unset = "TRUE"), "FALSE")

  if (is.na(grid_step) || grid_step <= 0) {
    stop("`AGETMP_GRID_STEP` must be positive for Figure 1F AD-TMP heatmap generation.", call. = FALSE)
  }

  message("No AD-TMP grid matrix found; generating it from source data.")
  clinical_raw <- temporalCPSA::ageTMP_load_discovery_clinical(data_dir)
  clinical <- data.frame(
    id = temporalCPSA::ageTMP_normalize_sample_ids(clinical_raw$id),
    age = suppressWarnings(as.numeric(clinical_raw$cDisc_age)),
    sex = clinical_raw$cDisc_Gender,
    dasc = clinical_raw$cDisc_clinical_status_at_collection_event,
    stringsAsFactors = FALSE
  )
  clinical_fit <- clinical[
    !is.na(clinical$age) &
      clinical$age <= fit_age_max &
      clinical$sex %in% c("Male", "Female"),
    ,
    drop = FALSE
  ]

  matrix_cache <- new.env(parent = emptyenv())
  metadata_cache <- new.env(parent = emptyenv())
  load_modality_matrix <- function(modality) {
    if (exists(modality, envir = matrix_cache, inherits = FALSE)) {
      return(get(modality, envir = matrix_cache, inherits = FALSE))
    }
    message("Loading ", modality, " feature matrix...")
    if (identical(modality, "phospho")) {
      raw <- temporalCPSA::ageTMP_load_molecular(data_dir = data_dir, modality = modality)
      split <- temporalCPSA::ageTMP_split_annotation_matrix(
        raw,
        annotation_cols = seq_len(min(9, ncol(raw))),
        row_id = "Site"
      )
      mat <- split$matrix
      phospho_gene <- split$annotation$ApprovedGeneSymbol
      names(phospho_gene) <- split$annotation$Site
      assign("phospho_gene", phospho_gene, envir = metadata_cache)
    } else {
      mat <- temporalCPSA::ageTMP_load_feature_matrix(
        data_dir = data_dir,
        modality = modality,
        collapse = TRUE,
        row_id = "ApprovedGeneSymbol"
      )
    }
    colnames(mat) <- temporalCPSA::ageTMP_normalize_sample_ids(colnames(mat))
    assign(modality, mat, envir = matrix_cache)
    mat
  }

  fit_tmp_matrix <- function(modality, prediction_metadata) {
    mat <- load_modality_matrix(modality)
    mat <- mat[, intersect(colnames(mat), clinical_fit$id), drop = FALSE]
    if (ncol(mat) < 5) {
      stop("Too few clinical samples match ", modality, " matrix.", call. = FALSE)
    }
    features <- rownames(mat)
    if (is.finite(feature_limit) && length(features) > feature_limit) {
      raw_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
      features <- names(sort(raw_sd, decreasing = TRUE))[seq_len(feature_limit)]
      message("Preselected top ", length(features), " raw-dynamic ", modality, " features for this run.")
    }

    message("Fitting ", modality, " AD-TMP grid with span ", span_by_modality[[modality]], "...")
    pred <- temporalCPSA::ageTMP_predict_tumor_trajectory_matrix(
      tumor_mat = mat,
      tumor_metadata = clinical_fit,
      prediction_metadata = prediction_metadata,
      features = features,
      tumor_sample_col = "id",
      tumor_age_col = "age",
      tumor_sex_col = "sex",
      prediction_sample_col = "id",
      prediction_age_col = "age",
      prediction_sex_col = "sex",
      center_age_range = c(0, center_age_max),
      fit_age_range = c(0, fit_age_max),
      span = span_by_modality[[modality]],
      adaptive_span = FALSE,
      n_cores = n_cores,
      progress = TRUE,
      prediction_scope = "all_samples",
      return_trajectory = FALSE
    )
    message("Finished estimating ", modality, " AD-TMP grid trajectories.")
    pred$matrix
  }

  grid_ages <- seq(0, fit_age_max, by = grid_step)
  grid_metadata <- do.call(rbind, lapply(c("Female", "Male"), function(sex) {
    data.frame(
      id = paste0(sex, "_age_", sprintf("%05.2f", grid_ages)),
      age = grid_ages,
      sex = sex,
      stringsAsFactors = FALSE
    )
  }))
  grid_tmp_matrices <- list(
    protein = fit_tmp_matrix("protein", grid_metadata),
    rna = fit_tmp_matrix("rna", grid_metadata),
    phospho = fit_tmp_matrix("phospho", grid_metadata)
  )
  grid_tmp_matrices <- lapply(grid_tmp_matrices, fill_edge_na_by_row)
  combined_grid <- temporalCPSA::ageTMP_combine_tmp_matrices(
    grid_tmp_matrices,
    proportion = dynamic_proportion,
    max_features = max_dynamic_features,
    common_samples = grid_metadata$id,
    scale_rows = TRUE
  )
  utils::write.table(
    data.frame(feature = rownames(combined_grid), combined_grid, check.names = FALSE),
    file = file.path(out_dir, "figure1f_ad_tmp_grid_tmp.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    grid_metadata,
    file = file.path(out_dir, "figure1f_ad_tmp_grid_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

generate_ad_tmp_sample_matrix <- function(data_dir, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fit_age_max <- env_num("AGETMP_FIT_AGE_MAX", 80)
  center_age_max <- env_num("AGETMP_CENTER_AGE_MAX", 80)
  span <- env_num("AGETMP_TRAJECTORY_SPAN", 1.5)
  span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_TRAJECTORY_SPAN", span),
    rna = env_num("AGETMP_RNA_TRAJECTORY_SPAN", span),
    phospho = env_num("AGETMP_PHOSPHO_TRAJECTORY_SPAN", span)
  )
  adaptive_span <- !identical(Sys.getenv("AGETMP_ADAPTIVE_TRAJECTORY_SPAN", unset = "TRUE"), "FALSE")
  min_span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_MIN_TRAJECTORY_SPAN", 0.5),
    rna = env_num("AGETMP_RNA_MIN_TRAJECTORY_SPAN", 1.0),
    phospho = env_num("AGETMP_PHOSPHO_MIN_TRAJECTORY_SPAN", 0.5)
  )
  max_span_by_modality <- c(
    protein = env_num("AGETMP_PROTEIN_MAX_TRAJECTORY_SPAN", 3.0),
    rna = env_num("AGETMP_RNA_MAX_TRAJECTORY_SPAN", 3.0),
    phospho = env_num("AGETMP_PHOSPHO_MAX_TRAJECTORY_SPAN", 3.0)
  )
  span_step <- env_num("AGETMP_TRAJECTORY_SPAN_STEP", 0.1)
  feature_limit <- env_num("AGETMP_FEATURE_LIMIT", Inf)
  dynamic_proportion <- env_num("AGETMP_DYNAMIC_PROPORTION", 0.5)
  max_dynamic_features <- env_num("AGETMP_MAX_DYNAMIC_FEATURES", 5000)
  figure1f_max_features_by_modality <- c(
    protein = env_num("AGETMP_FIGURE1F_PROTEIN_FEATURES", 5000),
    phospho = env_num("AGETMP_FIGURE1F_PHOSPHO_FEATURES", 3804),
    rna = env_num("AGETMP_FIGURE1F_RNA_FEATURES", 5000)
  )
  n_cores <- as.integer(env_num("AGETMP_N_CORES", 1))
  remove_dasc <- !identical(Sys.getenv("AGETMP_REMOVE_DASC", unset = "TRUE"), "FALSE")
  manuscript_style_rows <- !identical(Sys.getenv("AGETMP_FIGURE1F_MANUSCRIPT_STYLE_ROWS", unset = "TRUE"), "FALSE")
  balance_trajectory_sex <- !identical(Sys.getenv("AGETMP_BALANCE_TRAJECTORY_SEX", unset = "FALSE"), "FALSE")
  prediction_scope <- Sys.getenv("AGETMP_PREDICTION_SCOPE", "all_samples")
  legacy_figure1f_structure <- !identical(Sys.getenv("AGETMP_FIGURE1F_LEGACY_STRUCTURE", unset = "TRUE"), "FALSE")

  message("No study-sample AD-TMP matrix found; generating it from source data.")
  clinical_raw <- temporalCPSA::ageTMP_load_discovery_clinical(data_dir)
  clinical <- data.frame(
    id = temporalCPSA::ageTMP_normalize_sample_ids(clinical_raw$id),
    age = suppressWarnings(as.numeric(clinical_raw$cDisc_age)),
    sex = clinical_raw$cDisc_Gender,
    dasc = clinical_raw$cDisc_clinical_status_at_collection_event,
    stringsAsFactors = FALSE
  )
  clinical_fit <- clinical[
    !is.na(clinical$age) &
      clinical$age <= fit_age_max &
      clinical$sex %in% c("Male", "Female"),
    ,
    drop = FALSE
  ]
  clinical_prediction <- clinical_fit
  if (isTRUE(remove_dasc)) {
    clinical_prediction <- clinical_prediction[!clinical_prediction$dasc %in% "Deceased-due to disease", , drop = FALSE]
  }
  source_sample_reference_matrix <- Sys.getenv("AGETMP_SOURCE_SAMPLE_REFERENCE_MATRIX", "")
  if (nzchar(source_sample_reference_matrix) && file.exists(source_sample_reference_matrix)) {
    target_ids <- colnames(read_feature_matrix_file(source_sample_reference_matrix))
    target_ids <- target_ids[target_ids %in% clinical_prediction$id]
    clinical_prediction <- clinical_prediction[match(target_ids, clinical_prediction$id), , drop = FALSE]
    message(
      "Restricting source-derived sample predictions to ",
      nrow(clinical_prediction),
      " samples from `AGETMP_SOURCE_SAMPLE_REFERENCE_MATRIX`."
    )
  }

  matrix_cache <- new.env(parent = emptyenv())
  metadata_cache <- new.env(parent = emptyenv())
  load_modality_matrix <- function(modality) {
    if (exists(modality, envir = matrix_cache, inherits = FALSE)) {
      return(get(modality, envir = matrix_cache, inherits = FALSE))
    }
    message("Loading ", modality, " feature matrix...")
    if (identical(modality, "phospho")) {
      raw <- temporalCPSA::ageTMP_load_molecular(data_dir = data_dir, modality = modality)
      split <- temporalCPSA::ageTMP_split_annotation_matrix(
        raw,
        annotation_cols = seq_len(min(9, ncol(raw))),
        row_id = "Site"
      )
      mat <- split$matrix
      phospho_gene <- split$annotation$ApprovedGeneSymbol
      names(phospho_gene) <- split$annotation$Site
      assign("phospho_gene", phospho_gene, envir = metadata_cache)
    } else {
      mat <- temporalCPSA::ageTMP_load_feature_matrix(
        data_dir = data_dir,
        modality = modality,
        collapse = TRUE,
        row_id = "ApprovedGeneSymbol"
      )
    }
    colnames(mat) <- temporalCPSA::ageTMP_normalize_sample_ids(colnames(mat))
    assign(modality, mat, envir = matrix_cache)
    mat
  }

  modalities <- c("protein", "rna", "phospho")
  common_prediction_ids <- Reduce(
    intersect,
    c(list(clinical_prediction$id), lapply(modalities, function(modality) colnames(load_modality_matrix(modality))))
  )
  clinical_prediction <- clinical_prediction[match(common_prediction_ids, clinical_prediction$id), , drop = FALSE]
  message("Common plotted discovery prediction samples: ", nrow(clinical_prediction))

  read_matrix_cache <- function(path) {
    read_feature_matrix_file(path)
  }
  write_matrix_cache <- function(mat, path) {
    utils::write.table(
      data.frame(feature = rownames(mat), mat, check.names = FALSE),
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  read_fit_sample_reference <- function(modality) {
    env_name <- paste0("AGETMP_", toupper(modality), "_FIT_SAMPLE_REFERENCE")
    path <- Sys.getenv(env_name, "")
    if (!nzchar(path) || !file.exists(path)) {
      return(NULL)
    }
    ref <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if ("id" %in% names(ref)) {
      ids <- ref$id
    } else {
      ids <- ref[[1]]
    }
    ids <- temporalCPSA::ageTMP_normalize_sample_ids(ids)
    ids[nzchar(ids)]
  }
  write_fit_sample_cache <- function(ids, path) {
    utils::write.table(
      data.frame(id = ids),
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  fit_sample_cache_matches <- function(ids, path) {
    if (!file.exists(path)) {
      return(FALSE)
    }
    cached <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    cached_ids <- if ("id" %in% names(cached)) cached$id else cached[[1]]
    identical(as.character(cached_ids), as.character(ids))
  }
  collapse_phospho_tmp_by_gene <- function(mat) {
    if (!exists("phospho_gene", envir = metadata_cache, inherits = FALSE)) {
      stop("Phosphosite-to-gene map is missing.", call. = FALSE)
    }
    phospho_gene <- get("phospho_gene", envir = metadata_cache, inherits = FALSE)
    row_parts <- strsplit(rownames(mat), "::", fixed = TRUE)
    row_sex <- vapply(row_parts, function(x) {
      hit <- x[x %in% c("Female", "Male")]
      if (length(hit) > 0) hit[[1]] else NA_character_
    }, character(1))
    row_site <- vapply(row_parts, function(x) {
      hit <- x[!x %in% c("Female", "Male")]
      if (length(hit) > 0) hit[[1]] else x[[length(x)]]
    }, character(1))
    row_gene <- unname(phospho_gene[row_site])
    row_gene[is.na(row_gene) | !nzchar(row_gene)] <- row_site[is.na(row_gene) | !nzchar(row_gene)]
    row_group <- ifelse(
      !is.na(row_sex),
      paste(row_sex, row_gene, sep = "::"),
      row_gene
    )
    collapsed <- lapply(split(seq_len(nrow(mat)), row_group), function(i) {
      values <- colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
      values[is.nan(values)] <- NA_real_
      values
    })
    out <- do.call(rbind, collapsed)
    rownames(out) <- names(collapsed)
    colnames(out) <- colnames(mat)
    out
  }

  fill_legacy_prediction_na <- function(mat, metadata) {
    # Legacy phospho trajectories filled loess prediction gaps before
    # phosphosite-to-gene collapse, keeping neighboring age columns comparable.
    metadata <- metadata[match(colnames(mat), metadata$id), , drop = FALSE]
    sample_order <- order(metadata$age, na.last = TRUE)
    restore_order <- order(sample_order)
    filled <- t(apply(mat[, sample_order, drop = FALSE], 1, function(x) {
      if (!any(!is.na(x))) {
        return(x)
      }
      x[is.na(x)] <- stats::na.omit(x)[[1]]
      x
    }))
    colnames(filled) <- colnames(mat)[sample_order]
    filled[, restore_order, drop = FALSE]
  }

  fit_tmp_matrix <- function(modality, prediction_metadata) {
    cache_path <- file.path(
      out_dir,
      if (identical(modality, "phospho")) {
        paste0("figure1f_ad_tmp_", modality, "_site_full_tmp.tsv")
      } else {
        paste0("figure1f_ad_tmp_", modality, "_full_tmp.tsv")
      }
    )
    collapsed_cache_path <- file.path(out_dir, paste0("figure1f_ad_tmp_", modality, "_full_tmp.tsv"))
    cache_fit_samples_path <- file.path(out_dir, paste0("figure1f_ad_tmp_", modality, "_fit_samples.tsv"))
    mat <- load_modality_matrix(modality)
    fit_ids <- intersect(colnames(mat), clinical_fit$id)
    fit_reference_ids <- read_fit_sample_reference(modality)
    if (!is.null(fit_reference_ids)) {
      fit_ids <- fit_reference_ids[fit_reference_ids %in% fit_ids]
      message("Restricting ", modality, " trajectory fit to ", length(fit_ids), " reference samples.")
    } else {
      message("Using ", length(fit_ids), " available ", modality, " samples for trajectory fit.")
    }
    if (file.exists(cache_path) && fit_sample_cache_matches(fit_ids, cache_fit_samples_path)) {
      message("Loading cached ", modality, " study-sample AD-TMP matrix...")
      out <- read_matrix_cache(cache_path)
      if (identical(modality, "phospho")) {
        out <- fill_legacy_prediction_na(out, prediction_metadata)
        write_matrix_cache(out, cache_path)
        writeLines("phospho_legacy_prediction_na_fill_v1", file.path(out_dir, "figure1f_ad_tmp_phospho_legacy_na_fill.version"))
        load_modality_matrix("phospho")
        out <- collapse_phospho_tmp_by_gene(out)
        write_matrix_cache(out, collapsed_cache_path)
      }
      message("Finished loading cached ", modality, " study-sample AD-TMP matrix.")
      return(out)
    } else if (file.exists(cache_path)) {
      message("Cached ", modality, " AD-TMP matrix was fit on a different sample set; regenerating.")
    }
    mat <- mat[, fit_ids, drop = FALSE]
    if (ncol(mat) < 5) {
      stop("Too few clinical samples match ", modality, " matrix.", call. = FALSE)
    }
    features <- rownames(mat)
    if (is.finite(feature_limit) && length(features) > feature_limit) {
      raw_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
      features <- names(sort(raw_sd, decreasing = TRUE))[seq_len(feature_limit)]
      message("Preselected top ", length(features), " raw-dynamic ", modality, " features for this run.")
    }

    message(
      "Fitting ", modality, " study-sample AD-TMP matrix with ",
      if (isTRUE(adaptive_span)) {
        paste0("adaptive span ", min_span_by_modality[[modality]], "-", max_span_by_modality[[modality]])
      } else {
        paste0("span ", span_by_modality[[modality]])
      },
      "..."
    )
    pred <- temporalCPSA::ageTMP_predict_tumor_trajectory_matrix(
      tumor_mat = mat,
      tumor_metadata = clinical_fit,
      prediction_metadata = prediction_metadata,
      features = features,
      tumor_sample_col = "id",
      tumor_age_col = "age",
      tumor_sex_col = "sex",
      prediction_sample_col = "id",
      prediction_age_col = "age",
      prediction_sex_col = "sex",
      center_age_range = c(0, center_age_max),
      fit_age_range = c(0, fit_age_max),
      pre_scale = TRUE,
      span = span_by_modality[[modality]],
      adaptive_span = adaptive_span,
      min_span = min_span_by_modality[[modality]],
      max_span = max_span_by_modality[[modality]],
      span_step = span_step,
      n_cores = n_cores,
      progress = TRUE,
      prediction_scope = prediction_scope,
      return_trajectory = FALSE
    )
    out <- pred$matrix
    if (identical(modality, "phospho")) {
      out <- fill_legacy_prediction_na(out, prediction_metadata)
      write_matrix_cache(out, cache_path)
      writeLines("phospho_legacy_prediction_na_fill_v1", file.path(out_dir, "figure1f_ad_tmp_phospho_legacy_na_fill.version"))
      out <- collapse_phospho_tmp_by_gene(out)
      write_matrix_cache(out, collapsed_cache_path)
    } else {
      write_matrix_cache(out, cache_path)
    }
    write_fit_sample_cache(fit_ids, cache_fit_samples_path)
    message("Finished estimating ", modality, " study-sample AD-TMP trajectories.")
    out
  }

  row_trajectory_sex <- function(rows) {
    parts <- strsplit(rows, "::", fixed = TRUE)
    vapply(parts, function(x) {
      hit <- x[x %in% c("Female", "Male")]
      if (length(hit) > 0) hit[[1]] else NA_character_
    }, character(1))
  }

  tmp_matrices <- list(
    protein = fit_tmp_matrix("protein", clinical_prediction),
    rna = fit_tmp_matrix("rna", clinical_prediction),
    phospho = fit_tmp_matrix("phospho", clinical_prediction)
  )
  if (!isTRUE(legacy_figure1f_structure)) {
    tmp_matrices <- lapply(tmp_matrices, fill_edge_na_by_row)
  }

  round_even <- function(x) {
    out <- round(x)
    ifelse(out %% 2 == 0, out, out + 1 - 2 * (out %% 2))
  }
  figure1f_feature_name <- function(rows, modality) {
    rows <- as.character(rows)
    label <- if (identical(modality, "rna")) "trans" else modality
    parts <- strsplit(rows, "::", fixed = TRUE)
    vapply(parts, function(x) {
      sex <- x[x %in% c("Female", "Male")]
      feature <- x[!x %in% c("Female", "Male")]
      if (length(sex) == 0) sex <- NA_character_
      if (length(feature) == 0) feature <- x[length(x)]
      paste(sex[[1]], feature[[1]], label, sep = ",")
    }, character(1))
  }
  select_dynamic_rows <- function(mat, modality) {
    mat <- mat[, clinical_prediction$id, drop = FALSE]
    if (isTRUE(legacy_figure1f_structure)) {
      keep <- rowSums(!is.na(mat)) == ncol(mat)
      mat <- mat[keep, , drop = FALSE]
      row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
      row_sd[is.na(row_sd)] <- -Inf
      modality_max <- figure1f_max_features_by_modality[[modality]]
      n_keep <- min(modality_max, max_dynamic_features, max(1L, round_even(nrow(mat) / 2)))
      mat <- mat[order(row_sd, decreasing = TRUE)[seq_len(n_keep)], , drop = FALSE]
      mat <- t(scale(t(mat)))
      mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
      rownames(mat) <- make.unique(figure1f_feature_name(rownames(mat), modality))
      return(mat)
    }
    keep <- rowSums(!is.na(mat)) >= 2
    mat <- mat[keep, , drop = FALSE]
    row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
    row_sd[is.na(row_sd)] <- -Inf
    modality_max <- figure1f_max_features_by_modality[[modality]]
    n_keep <- min(modality_max, max_dynamic_features, max(1L, round_even(nrow(mat) / 2)))
    mat <- mat[order(row_sd, decreasing = TRUE)[seq_len(n_keep)], , drop = FALSE]
    mat <- t(scale(t(mat)))
    mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
    rownames(mat) <- make.unique(figure1f_feature_name(rownames(mat), modality))
    mat
  }

  if (isTRUE(manuscript_style_rows)) {
    selected <- lapply(c("rna", "protein", "phospho"), function(modality) {
      select_dynamic_rows(tmp_matrices[[modality]], modality)
    })
    names(selected) <- c("rna", "protein", "phospho")
    if (isTRUE(legacy_figure1f_structure)) {
      common_ids <- Reduce(intersect, lapply(selected, colnames))
      selected <- lapply(selected, function(mat) mat[, match(common_ids, colnames(mat)), drop = FALSE])
    }
    combined <- do.call(rbind, selected)
    if (isTRUE(legacy_figure1f_structure)) {
      combined <- t(scale(t(combined)))
      combined <- combined[rowSums(is.na(combined)) == 0, , drop = FALSE]
      clinical_prediction <- clinical_prediction[match(colnames(combined), clinical_prediction$id), , drop = FALSE]
    }
  } else if (isTRUE(balance_trajectory_sex)) {
    selected <- lapply(names(tmp_matrices), function(modality) {
      mat <- tmp_matrices[[modality]][, clinical_prediction$id, drop = FALSE]
      sex <- row_trajectory_sex(rownames(mat))
      per_sex_max <- max(1L, floor(max_dynamic_features / 2))
      by_sex <- lapply(c("Female", "Male"), function(sex_level) {
        sub <- mat[sex %in% sex_level, , drop = FALSE]
        if (nrow(sub) == 0) {
          return(sub)
        }
        temporalCPSA::ageTMP_select_dynamic_tmp_features(
          sub,
          proportion = dynamic_proportion,
          max_features = per_sex_max,
          scale_rows = TRUE
        )
      })
      mat <- do.call(rbind, by_sex)
      dimnames(mat) <- list(paste(modality, make.unique(as.character(rownames(mat))), sep = "::"), colnames(mat))
      mat
    })
    names(selected) <- names(tmp_matrices)
    combined <- do.call(rbind, selected)
  } else {
    combined <- temporalCPSA::ageTMP_combine_tmp_matrices(
      tmp_matrices,
      proportion = dynamic_proportion,
      max_features = max_dynamic_features,
      common_samples = clinical_prediction$id,
      scale_rows = TRUE
    )
  }

  combined <- combined[, intersect(clinical_prediction$id, colnames(combined)), drop = FALSE]
  utils::write.table(
    data.frame(feature = rownames(combined), combined, check.names = FALSE),
    file = file.path(out_dir, "figure1f_ad_tmp_combined_tmp.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    clinical_prediction,
    file = file.path(out_dir, "figure1f_ad_tmp_sample_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

figure1f_display_mode <- Sys.getenv(
  "AGETMP_FIGURE1F_DISPLAY_MODE",
  "samples"
)
if (!figure1f_display_mode %in% c("samples", "grid")) {
  stop("`AGETMP_FIGURE1F_DISPLAY_MODE` must be `samples` or `grid`.", call. = FALSE)
}
legacy_clustme_rdata <- Sys.getenv("AGETMP_LEGACY_CLUSTME_RDATA", "")
legacy_clustme_matrix <- Sys.getenv(
  "AGETMP_LEGACY_CLUSTME_MATRIX",
  file.path(script_dir, "legacy_inputs", "figure1f_ad_tmp_legacy_clustme.tsv")
)
legacy_clustme_output_matrix <- file.path(out_dir, "figure1f_ad_tmp_legacy_clustme.tsv")
if (nzchar(legacy_clustme_rdata) && file.exists(legacy_clustme_rdata) && !file.exists(legacy_clustme_matrix)) {
  legacy_env <- new.env(parent = emptyenv())
  load(legacy_clustme_rdata, envir = legacy_env)
  if (!exists("comb.age.col.clusters", envir = legacy_env, inherits = FALSE)) {
    stop("Legacy clustme RData does not contain `comb.age.col.clusters`.", call. = FALSE)
  }
  legacy_clusters <- get("comb.age.col.clusters", envir = legacy_env)
  if (length(legacy_clusters) < 3 || is.null(legacy_clusters[[3]]$clustme)) {
    stop("Legacy clustme RData does not contain `comb.age.col.clusters[[3]]$clustme`.", call. = FALSE)
  }
  legacy_clustme <- legacy_clusters[[3]]$clustme
  utils::write.table(
    data.frame(feature = rownames(legacy_clustme), legacy_clustme, check.names = FALSE),
    file = legacy_clustme_matrix,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
if (file.exists(legacy_clustme_matrix) && !file.exists(legacy_clustme_output_matrix)) {
  utils::write.table(
    utils::read.delim(legacy_clustme_matrix, check.names = FALSE),
    file = legacy_clustme_output_matrix,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
tmp_matrix_candidates <- if (identical(figure1f_display_mode, "grid")) {
  c(
    Sys.getenv("AGETMP_TMP_MATRIX", ""),
    file.path(out_dir, "figure1f_ad_tmp_grid_tmp.tsv"),
    file.path(out_dir, "figure1f_ad_tmp_combined_tmp.tsv"),
    file.path(out_dir, "figure1f_ad_tmp_smoke_combined_tmp.tsv")
  )
} else {
  if (identical(figure1f_mode, "manuscript")) {
    c(
      Sys.getenv("AGETMP_TMP_MATRIX", ""),
      legacy_clustme_matrix,
      legacy_clustme_output_matrix,
      file.path(out_dir, "figure1f_legacy_clustme", "figure1f_ad_tmp_legacy_clustme.tsv")
    )
  } else {
    c(
      Sys.getenv("AGETMP_TMP_MATRIX", ""),
      file.path(out_dir, "figure1f_ad_tmp_combined_tmp.tsv"),
      file.path(out_dir, "figure1f_ad_tmp_smoke_combined_tmp.tsv")
    )
  }
}
tmp_matrix_path <- tmp_matrix_candidates[nzchar(tmp_matrix_candidates) & file.exists(tmp_matrix_candidates)][1]
source_fit_reference_paths <- vapply(
  c("protein", "rna", "phospho"),
  function(modality) Sys.getenv(paste0("AGETMP_", toupper(modality), "_FIT_SAMPLE_REFERENCE"), ""),
  character(1)
)
source_fit_cache_paths <- file.path(
  out_dir,
  paste0("figure1f_ad_tmp_", c("protein", "rna", "phospho"), "_fit_samples.tsv")
)
source_fit_references_active <- any(nzchar(source_fit_reference_paths) & file.exists(source_fit_reference_paths))
source_fit_cache_complete <- all(file.exists(source_fit_cache_paths))
if (
  identical(figure1f_mode, "source") &&
    !identical(figure1f_display_mode, "grid") &&
    !is.na(tmp_matrix_path) &&
    identical(normalizePath(tmp_matrix_path, mustWork = FALSE), normalizePath(file.path(out_dir, "figure1f_ad_tmp_combined_tmp.tsv"), mustWork = FALSE)) &&
    (
      !file.exists(file.path(out_dir, "figure1f_ad_tmp_phospho_site_full_tmp.tsv")) ||
        !file.exists(file.path(out_dir, "figure1f_ad_tmp_phospho_legacy_na_fill.version")) ||
        (source_fit_references_active && !source_fit_cache_complete)
    )
) {
  message("Existing source-derived Figure 1F matrix predates current phospho trajectory handling; regenerating.")
  tmp_matrix_path <- NA_character_
}
if (is.na(tmp_matrix_path) || !nzchar(tmp_matrix_path)) {
  if (identical(figure1f_display_mode, "grid")) {
    generate_ad_tmp_grid(data_dir = data_dir, out_dir = out_dir)
    tmp_matrix_path <- file.path(out_dir, "figure1f_ad_tmp_grid_tmp.tsv")
  } else if (identical(figure1f_mode, "manuscript")) {
    stop(
      "No archived Figure 1F AD-TMP matrix was found. ",
      "Place `figure1f_ad_tmp_legacy_clustme.tsv` in `legacy_inputs/`, ",
      "set `AGETMP_LEGACY_CLUSTME_MATRIX`, or run with ",
      "`AGETMP_FIGURE1F_MODE=source` for a source-derived analogous plot.",
      call. = FALSE
    )
  } else {
    generate_ad_tmp_sample_matrix(data_dir = data_dir, out_dir = out_dir)
    tmp_matrix_path <- file.path(out_dir, "figure1f_ad_tmp_combined_tmp.tsv")
  }
  if (!file.exists(tmp_matrix_path)) {
    stop("AD-TMP generation did not create the expected matrix file.", call. = FALSE)
  }
}

read_tmp_matrix <- function(path) {
  if (
    identical(basename(path), "figure1f_ad_tmp_legacy_clustme.tsv") &&
      exists("ageTMP_load_figure1f_ad_tmp_matrix", where = asNamespace("temporalCPSA"), inherits = FALSE)
  ) {
    mat <- temporalCPSA::ageTMP_load_figure1f_ad_tmp_matrix(path)
  } else if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    mat <- readRDS(path)
  } else {
    mat <- utils::read.delim(path, check.names = FALSE)
    rownames(mat) <- mat[[1]]
    mat <- mat[, -1, drop = FALSE]
  }
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

write_tmp_matrix <- function(mat, path) {
  utils::write.table(
    data.frame(feature = rownames(mat), mat, check.names = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

clean_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

fill_edge_na <- function(x) {
  ok <- which(!is.na(x))
  if (length(ok) == 0) {
    return(x)
  }
  if (length(ok) == 1) {
    return(rep(x[ok], length(x)))
  }
  stats::approx(ok, x[ok], xout = seq_along(x), rule = 2)$y
}

matched_optional <- function(data, column, ids, id_column = "id") {
  if (!column %in% names(data)) {
    return(rep(NA, length(ids)))
  }
  data[[column]][match(ids, data[[id_column]])]
}

age_class <- as.data.frame(temporalCPSA::ageTMP_load_supplement(
  data_dir = data_dir,
  table = "STable1",
  sheet = "cDisc_AgeClass"
))
clinical <- as.data.frame(temporalCPSA::ageTMP_load_supplement(
  data_dir = data_dir,
  table = "STable1",
  sheet = "ClinicalTable"
))
validation <- as.data.frame(temporalCPSA::ageTMP_load_supplement(
  data_dir = data_dir,
  table = "STable1",
  sheet = "Val_ClinicalTable"
))
reference_clinical <- as.data.frame(temporalCPSA::ageTMP_load_supplement(
  data_dir = data_dir,
  table = "STable1",
  sheet = "Ref_ClinicalTable"
))

grid_metadata_file <- file.path(out_dir, "figure1f_ad_tmp_grid_metadata.tsv")
grid_mode <- FALSE

tmp_matrix <- read_tmp_matrix(tmp_matrix_path)
colnames(tmp_matrix) <- temporalCPSA::ageTMP_normalize_sample_ids(colnames(tmp_matrix))
grid_mode <- file.exists(grid_metadata_file) &&
  all(colnames(tmp_matrix) %in% temporalCPSA::ageTMP_normalize_sample_ids(utils::read.delim(grid_metadata_file, check.names = FALSE)$id))

if (isTRUE(grid_mode)) {
  grid_metadata <- utils::read.delim(grid_metadata_file, check.names = FALSE)
  grid_metadata$id <- temporalCPSA::ageTMP_normalize_sample_ids(grid_metadata$id)
  grid_metadata$age <- clean_numeric(grid_metadata$age)
  grid_metadata$sex <- as.character(grid_metadata$sex)
  ages <- sort(unique(grid_metadata$age))
  by_sex <- lapply(c("Female", "Male"), function(sex) {
    ids <- grid_metadata$id[grid_metadata$sex == sex][order(grid_metadata$age[grid_metadata$sex == sex])]
    mat <- tmp_matrix[, ids, drop = FALSE]
    colnames(mat) <- paste0("age_", sprintf("%05.2f", ages))
    rownames(mat) <- paste(sex, rownames(mat), sep = "::")
    mat
  })
  tmp_matrix <- do.call(rbind, by_sex)
  sample_meta <- data.frame(
    id = paste0("age_", sprintf("%05.2f", ages)),
    age = ages,
    age_class = factor(
      as.character(cut(ages, breaks = c(0, 15, 26, 40, 62, 80), include.lowest = TRUE, labels = c("PED", "ADO", "YA", "ADULT", "SEN"))),
      levels = c("PED", "ADO", "YA", "ADULT", "SEN"),
      ordered = TRUE
    ),
    published_sil = NA_real_,
    sex = NA_character_,
    H33A_mut = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  tmp_matrix <- tmp_matrix[, colnames(tmp_matrix) %in% temporalCPSA::ageTMP_normalize_sample_ids(age_class$id), drop = FALSE]
  if (ncol(tmp_matrix) < 5) {
    stop("Fewer than five age-class samples are present in the TMP matrix.", call. = FALSE)
  }
  sample_meta <- data.frame(
    id = temporalCPSA::ageTMP_normalize_sample_ids(age_class$id),
    age = clean_numeric(clinical$cDisc_age[match(age_class$id, clinical$id)]),
    age_class = factor(age_class$age.class, levels = c("PED", "ADO", "YA", "ADULT", "SEN"), ordered = TRUE),
    published_sil = clean_numeric(age_class$sil_width),
    sex = clinical$cDisc_Gender[match(age_class$id, clinical$id)],
    H33A_mut = suppressWarnings(as.numeric(matched_optional(clinical, "cDisc_H33A_mut", age_class$id))),
    stringsAsFactors = FALSE
  )
  sample_meta <- sample_meta[match(colnames(tmp_matrix), sample_meta$id), , drop = FALSE]
  sample_meta <- sample_meta[order(sample_meta$age), , drop = FALSE]
  tmp_matrix <- tmp_matrix[, sample_meta$id, drop = FALSE]
}

heatmap_feature_limit <- env_num("AGETMP_HEATMAP_FEATURES", Inf)
preserve_figure1f_row_order <- !identical(Sys.getenv("AGETMP_FIGURE1F_PRESERVE_ROW_ORDER", unset = "TRUE"), "FALSE")
if (isTRUE(preserve_figure1f_row_order) && !is.finite(heatmap_feature_limit)) {
  plotme <- tmp_matrix
} else {
  row_sd <- apply(tmp_matrix, 1, stats::sd, na.rm = TRUE)
  if (is.finite(heatmap_feature_limit)) {
  keep_n <- min(length(row_sd), as.integer(heatmap_feature_limit))
  keep_features <- names(sort(row_sd, decreasing = TRUE))[seq_len(keep_n)]
  plotme <- tmp_matrix[keep_features, , drop = FALSE]
  } else {
  plotme <- tmp_matrix[order(row_sd, decreasing = TRUE), , drop = FALSE]
  }
  plotme <- t(scale(t(plotme)))
}
dtt_limit <- as.numeric(Sys.getenv("AGETMP_HEATMAP_LIMIT", "2"))
plotme[plotme > dtt_limit] <- dtt_limit
plotme[plotme < -dtt_limit] <- -dtt_limit

match_plot_metadata <- function(plot_matrix, metadata) {
  if (!"id" %in% names(metadata)) {
    stop("Plot metadata must contain an `id` column.", call. = FALSE)
  }
  if (anyDuplicated(colnames(plot_matrix)) > 0) {
    stop("Plot matrix contains duplicated sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(metadata$id) > 0) {
    stop("Plot metadata contains duplicated sample IDs.", call. = FALSE)
  }

  keep <- colnames(plot_matrix) %in% metadata$id
  dropped <- setdiff(colnames(plot_matrix), metadata$id)
  if (length(dropped) > 0) {
    warning(
      "Dropping ", length(dropped),
      " plot matrix column(s) without matching metadata: ",
      paste(head(dropped, 5), collapse = ", "),
      if (length(dropped) > 5) ", ..." else "",
      call. = FALSE
    )
  }
  plot_matrix <- plot_matrix[, keep, drop = FALSE]
  metadata <- metadata[match(colnames(plot_matrix), metadata$id), , drop = FALSE]
  if (!identical(colnames(plot_matrix), metadata$id)) {
    stop("Failed to align plot matrix columns with sample metadata rows.", call. = FALSE)
  }
  list(matrix = plot_matrix, metadata = metadata)
}

matched_plot <- match_plot_metadata(plotme, sample_meta)
plotme <- matched_plot$matrix
sample_meta <- matched_plot$metadata

default_survival_source <- if (identical(figure1f_mode, "manuscript")) {
  "figure1f_legacy"
} else {
  "ref"
}
survival_source <- Sys.getenv("AGETMP_SURVIVAL_SOURCE", default_survival_source)
legacy_survival_annotation <- Sys.getenv(
  "AGETMP_FIGURE1F_SURVIVAL_ANNOTATION",
  file.path(script_dir, "legacy_inputs", "figure1f_survival_annotation_data.tsv")
)
sample_meta$m_os_vali <- NA_real_
sample_meta$f_os_vali <- NA_real_
if (survival_source == "figure1f_legacy") {
  if (!file.exists(legacy_survival_annotation)) {
    stop(
      "No archived Figure 1F survival-days annotation file was found. ",
      "Place `figure1f_survival_annotation_data.tsv` in `legacy_inputs/`, ",
      "set `AGETMP_FIGURE1F_SURVIVAL_ANNOTATION`, or set ",
      "`AGETMP_SURVIVAL_SOURCE=ref` for a source-derived survival annotation.",
      call. = FALSE
    )
  }
  figure1f_survival <- temporalCPSA::ageTMP_load_figure1f_survival_annotation(legacy_survival_annotation)
  survival_for_loess <- data.frame(
    id = figure1f_survival$id,
    age_numeric = clean_numeric(figure1f_survival$age),
    days = clean_numeric(figure1f_survival$days),
    os_status_numeric = clean_numeric(figure1f_survival$os.status),
    sex = as.character(figure1f_survival$Gender),
    stringsAsFactors = FALSE
  )
} else if (survival_source == "ref") {
  survival_for_loess <- data.frame(
    id = reference_clinical$id,
    age_numeric = clean_numeric(reference_clinical$age),
    days = clean_numeric(reference_clinical$days),
    os_status_numeric = clean_numeric(reference_clinical$os.status),
    sex = as.character(reference_clinical$Gender),
    stringsAsFactors = FALSE
  )
} else if (survival_source == "validation") {
  survival_for_loess <- data.frame(
    id = validation$id,
    age_numeric = clean_numeric(validation$age),
    days = clean_numeric(validation$survival_days),
    os_status_numeric = clean_numeric(validation$os_status),
    sex = as.character(validation$Gender),
    stringsAsFactors = FALSE
  )
} else if (survival_source == "cdisc") {
  survival_for_loess <- data.frame(
    id = clinical$id,
    age_numeric = clean_numeric(clinical$cDisc_age),
    days = clean_numeric(clinical$cDisc_os),
    os_status_numeric = clean_numeric(clinical$cDisc_os_status),
    sex = as.character(clinical$cDisc_Gender),
    stringsAsFactors = FALSE
  )
} else {
  stop("`AGETMP_SURVIVAL_SOURCE` must be `figure1f_legacy`, `ref`, `validation`, or `cdisc`.", call. = FALSE)
}

loess_data <- survival_for_loess[
  survival_for_loess$os_status_numeric %in% 1 &
    !is.na(survival_for_loess$days) &
    !is.na(survival_for_loess$age_numeric),
  ,
  drop = FALSE
]
if (nrow(loess_data) < 5) {
  stop("Too few validation event samples to fit the LOESS survival annotation.", call. = FALSE)
}
loess_span <- as.numeric(Sys.getenv("AGETMP_SURVIVAL_LOESS_SPAN", "0.30"))
loess_vali_os <- stats::loess(days ~ age_numeric, data = loess_data, span = loess_span)
sample_meta$os_vali_raw <- as.numeric(predict(loess_vali_os, newdata = data.frame(age_numeric = sample_meta$age)))
sample_meta$os_vali <- sample_meta$os_vali_raw

female_data <- loess_data[loess_data$sex == "Female", , drop = FALSE]
male_data <- loess_data[loess_data$sex == "Male", , drop = FALSE]
if (nrow(female_data) >= 5) {
  female_fit <- stats::loess(days ~ age_numeric, data = female_data, span = loess_span)
  sample_meta$f_os_vali <- as.numeric(predict(female_fit, newdata = data.frame(age_numeric = sample_meta$age)))
}
if (nrow(male_data) >= 5) {
  male_fit <- stats::loess(days ~ age_numeric, data = male_data, span = loess_span)
  sample_meta$m_os_vali <- as.numeric(predict(male_fit, newdata = data.frame(age_numeric = sample_meta$age)))
}
survival_edge_fill <- !identical(Sys.getenv("AGETMP_SURVIVAL_EDGE_FILL", "TRUE"), "FALSE")
if (survival_edge_fill && any(is.na(sample_meta$os_vali))) {
  sample_meta$os_vali <- fill_edge_na(sample_meta$os_vali)
}
survival_clamp_zero <- !identical(Sys.getenv("AGETMP_SURVIVAL_CLAMP_ZERO", "TRUE"), "FALSE")
if (survival_clamp_zero) {
  sample_meta$os_vali <- pmax(0, sample_meta$os_vali)
}

sample_meta$sil_width <- NA_real_
if (requireNamespace("cluster", quietly = TRUE) && all(table(sample_meta$age_class) >= 2)) {
  sil <- cluster::silhouette(as.numeric(sample_meta$age_class), stats::dist(t(plotme)))
  sample_meta$sil_width <- sil[, "sil_width"]
}

age_class_col <- c(
  PED = "#c7e9c0",
  ADO = "#a1d99b",
  YA = "#74c476",
  ADULT = "#238b45",
  SEN = "#005a32"
)
sex_col <- c(Female = "#c65a9b", Male = "#3c78b5")
parse_row_tokens <- function(rows) {
  lapply(rows, function(row) {
    if (grepl("::", row, fixed = TRUE)) {
      strsplit(row, "::", fixed = TRUE)[[1]]
    } else if (grepl(",", row, fixed = TRUE)) {
      strsplit(row, ",", fixed = TRUE)[[1]]
    } else {
      row
    }
  })
}
row_parts <- parse_row_tokens(rownames(plotme))
modality <- vapply(row_parts, function(x) {
  hit <- x[x %in% c("protein", "rna", "trans", "phospho")]
  if (length(hit) > 0) {
    if (identical(hit[[1]], "trans")) "rna" else hit[[1]]
  } else {
    NA_character_
  }
}, character(1))
row_sex <- vapply(row_parts, function(x) {
  hit <- x[x %in% c("Female", "Male")]
  if (length(hit) > 0) hit[[1]] else NA_character_
}, character(1))
row_sex[!row_sex %in% c("Female", "Male")] <- NA_character_
modality_col <- c(rna = "#a020f0", protein = "#ffa500", phospho = "#00ff00")
sex_col <- c(Female = "#cc0000", Male = "#0000cc")
row_split <- factor(
  paste(modality, row_sex, sep = " / "),
  levels = c(
    "rna / Male",
    "rna / Female",
    "protein / Male",
    "protein / Female",
    "phospho / Male",
    "phospho / Female"
  )
)

top_anno <- HeatmapAnnotation(
  survival = anno_lines(
    sample_meta[, "os_vali", drop = FALSE],
    gp = gpar(col = "#35dce8", lwd = 1.4),
    smooth = FALSE,
    add_points = TRUE,
    ylim = c(0, max(1200, sample_meta$os_vali, na.rm = TRUE)),
    pt_gp = gpar(col = "#35dce8", lwd = 1),
    height = unit(2.3, "cm"),
    pch = 1,
    axis_param = list(
      side = "left",
      at = c(0, 500, 1000),
      gp = gpar(fontsize = 6)
    )
  ),
  age = anno_barplot(
    sample_meta$age,
    ylim = c(0, 80),
    gp = gpar(col = NA, fill = age_class_col[as.character(sample_meta$age_class)]),
    height = unit(0.75, "cm"),
    border = FALSE,
    axis_param = list(
      side = "left",
      at = c(0, 40, 80),
      gp = gpar(fontsize = 6)
    )
  ),
  sil = anno_barplot(
    sample_meta$sil_width,
    ylim = c(-0.5, 0.8),
    gp = gpar(col = NA, fill = "gray60"),
    height = unit(0.75, "cm"),
    border = FALSE,
    axis_param = list(
      side = "left",
      at = c(-0.4, 0, 0.4, 0.8),
      gp = gpar(fontsize = 6)
    )
  ),
  na_col = "white",
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8)
)

right_anno <- rowAnnotation(
  `Data type` = modality,
  sex = row_sex,
  col = list(`Data type` = modality_col, sex = sex_col),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm")
)

ht <- Heatmap(
  plotme,
  name = "AD-TMP",
  col = colorRamp2(c(-dtt_limit, 0, dtt_limit), c("blue", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  column_split = sample_meta$age_class,
  row_split = row_split,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  row_gap = unit(0, "mm"),
  use_raster = TRUE,
  width = unit(3.6, "in"),
  height = unit(5.6, "in"),
  top_annotation = top_anno,
  right_annotation = right_anno,
  row_title = NULL,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 14),
  heatmap_legend_param = list(
    at = c(-dtt_limit, -dtt_limit / 2, 0, dtt_limit / 2, dtt_limit),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

png_file <- file.path(out_dir, paste0(output_prefix, ".png"))
pdf_file <- file.path(out_dir, paste0(output_prefix, ".pdf"))
png_type <- if (capabilities("cairo")) "cairo" else getOption("bitmapType", "quartz")
png(png_file, width = 1450, height = 2100, res = 260, type = png_type)
ht_drawn <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
pdf(pdf_file, width = 5.6, height = 8.1)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

drawn_row_order <- unlist(row_order(ht_drawn), use.names = FALSE)
plotme_drawn <- plotme[drawn_row_order, , drop = FALSE]
modality_drawn <- modality[drawn_row_order]
row_sex_drawn <- row_sex[drawn_row_order]
row_split_drawn <- row_split[drawn_row_order]
utils::write.table(
  data.frame(
    row_index = seq_len(nrow(plotme_drawn)),
    feature = rownames(plotme_drawn),
    modality = modality_drawn,
    sex = row_sex_drawn,
    row_split = as.character(row_split_drawn),
    stringsAsFactors = FALSE
  ),
  file = file.path(out_dir, paste0(output_prefix, "_row_annotations.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  data.frame(feature = rownames(plotme_drawn), plotme_drawn, check.names = FALSE),
  file = file.path(out_dir, paste0(output_prefix, "_ordered_matrix.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

utils::write.table(
  sample_meta,
  file = file.path(out_dir, paste0(output_prefix, "_sample_annotations.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  sample_meta[, c("id", "age", "age_class", "os_vali_raw", "os_vali", "f_os_vali", "m_os_vali"), drop = FALSE],
  file = file.path(out_dir, paste0(output_prefix, "_survival_grid.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
summary_df <- data.frame(
  metric = c(
    "tmp_matrix",
    "figure1f_mode",
    "display_mode",
    "n_samples",
    "n_tmp_features_plotted",
    "n_survival_event_samples_for_loess",
    "loess_span",
    "grid_mode",
    "survival_source",
    "survival_edge_fill",
    "survival_clamp_zero",
    "survival_raw_na_count",
    "heatmap_feature_limit",
    "row_order_method",
    "tmp_note"
  ),
  value = c(
    tmp_matrix_path,
    figure1f_mode,
    figure1f_display_mode,
    ncol(plotme),
    nrow(plotme),
    nrow(loess_data),
    loess_span,
    grid_mode,
    survival_source,
    survival_edge_fill,
    survival_clamp_zero,
    sum(is.na(sample_meta$os_vali_raw)),
    heatmap_feature_limit,
    "source-derived AD-TMP rows use manuscript-style feature counts and ComplexHeatmap clustering; silhouette is recomputed from this plotted matrix",
    if (identical(figure1f_mode, "manuscript")) {
      "uses archived manuscript clustme matrix for exact Figure 1F visualization reproduction"
    } else if (grepl("smoke", tmp_matrix_path)) {
      "preview uses source-derived smoke TMP matrix; run full source TMP clustering for final reproduction"
    } else {
      "uses source-derived TMP matrix for an analogous Figure 1F reanalysis"
    }
  )
)
utils::write.table(
  summary_df,
  file = file.path(out_dir, paste0(output_prefix, "_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Figure 1F AD-TMP heatmap written:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat("Mode:", figure1f_mode, "\n")
cat("TMP matrix:", tmp_matrix_path, "\n")
cat("Samples:", ncol(plotme), "\n")
cat("Survival source:", survival_source, "\n")
cat("Survival event samples for LOESS:", nrow(loess_data), "\n")
