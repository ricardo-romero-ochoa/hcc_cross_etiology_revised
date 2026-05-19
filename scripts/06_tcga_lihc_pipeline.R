source("R/_shared.R")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)
  library(survival)
  library(broom)
})

ensure_dirs()
panels <- read_module_panels()
panels <- lapply(panels, toupper)

# -----------------------------
# Helpers specific to TCGA/GDC
# -----------------------------
make_coldata_tibble <- function(se) {
  meta0 <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)

  # GDCprepare/TCGAbiolinks objects may already contain a sample_id column.
  # rownames_to_column("sample_id") fails in that case, so preserve any existing
  # field under a different name and use the SE column names as the canonical ID.
  if ("sample_id" %in% names(meta0)) {
    names(meta0)[names(meta0) == "sample_id"] <- "sample_id_coldata"
  }
  names(meta0) <- make.unique(names(meta0), sep = "_")

  meta <- tibble::rownames_to_column(meta0, var = "sample_id") |>
    tibble::as_tibble()

  # Some SE objects use simple rownames but expression colnames carry the TCGA barcode.
  # If possible, force metadata IDs to match the assay colnames exactly.
  if (nrow(meta) == ncol(se)) {
    meta$sample_id <- colnames(se)
  }
  meta
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

infer_tcga_tissue <- function(meta) {
  # Prefer the TCGA barcode sample code because it is stable across TCGAbiolinks versions:
  # 01 = primary tumor; 11 = solid tissue normal. Fall back to text metadata.
  sample_code <- suppressWarnings(substr(meta$sample_id, 14, 15))
  tissue_by_code <- dplyr::case_when(
    sample_code %in% c("01", "02", "03", "05", "06", "07", "08", "09") ~ "tumor",
    sample_code %in% c("10", "11", "12", "13", "14") ~ "non_tumor",
    TRUE ~ NA_character_
  )

  text_cols <- intersect(
    c("sample_type", "definition", "shortLetterCode", "tissue_type", "sample", "sample_id_coldata"),
    names(meta)
  )
  text <- rep("", nrow(meta))
  if (length(text_cols) > 0) {
    text <- apply(meta[, text_cols, drop = FALSE], 1, function(x) {
      paste(as.character(x[!is.na(x)]), collapse = " | ")
    })
  }

  tissue_by_text <- dplyr::case_when(
    stringr::str_detect(text, stringr::regex("solid tissue normal|normal tissue|adjacent|non[- ]?tumou?r|non[- ]?cancer", TRUE)) ~ "non_tumor",
    stringr::str_detect(text, stringr::regex("primary tumor|tumou?r|carcinoma|cancer", TRUE)) ~ "tumor",
    TRUE ~ NA_character_
  )

  dplyr::coalesce(tissue_by_code, tissue_by_text)
}

collapse_counts_to_symbols <- function(expr_counts, rowdata) {
  rd <- as.data.frame(rowdata, stringsAsFactors = FALSE)
  sym_col <- first_existing_col(rd, c("gene_name", "external_gene_name", "gene_symbol", "hgnc_symbol"))
  if (is.na(sym_col)) {
    stop("Could not find a gene-symbol column in TCGA rowData. Columns found: ", paste(names(rd), collapse = ", "))
  }

  sym <- toupper(as.character(rd[[sym_col]]))
  ok <- !is.na(sym) & sym != "" & !stringr::str_detect(sym, "^ENSG")
  expr_counts <- expr_counts[ok, , drop = FALSE]
  sym <- sym[ok]
  storage.mode(expr_counts) <- "numeric"
  rowsum(expr_counts, group = sym, reorder = FALSE)
}

clean_clinical_for_merge <- function(clin) {
  clin <- tibble::as_tibble(clin)
  if (!"submitter_id" %in% names(clin)) {
    stop("TCGA clinical table lacks submitter_id; cannot merge clinical metadata.")
  }
  clin |>
    dplyr::mutate(patient = .data$submitter_id) |>
    dplyr::distinct(.data$patient, .keep_all = TRUE)
}

fit_lm_associations <- function(scores2) {
  candidates <- c("ajcc_pathologic_stage", "tumor_grade", "grade", "gender", "sex", "age_at_index")
  score_names <- c("ProlifHubScore", "HepLossScore", "HCCStateScore")
  assoc <- list()

  for (v in intersect(candidates, names(scores2))) {
    for (s in score_names) {
      dat <- scores2 |>
        dplyr::filter(.data$tissue == "tumor", !is.na(.data[[v]]), !is.na(.data[[s]]))
      if (nrow(dat) < 10) next
      if (!is.numeric(dat[[v]])) {
        dat[[v]] <- as.factor(dat[[v]])
        if (nlevels(dat[[v]]) < 2) next
      }
      fit <- tryCatch(stats::lm(stats::reformulate(v, response = s), data = dat), error = function(e) NULL)
      if (!is.null(fit)) {
        assoc[[paste(v, s, sep = "__")]] <- broom::tidy(fit) |>
          dplyr::mutate(variable = v, score = s, n = nrow(dat))
      }
    }
  }

  if (length(assoc) == 0) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(assoc)
}

fit_survival_models <- function(scores2) {
  coalesce_cols <- function(df, cols) {
    cols <- cols[cols %in% names(df)]
    if (length(cols) == 0) return(rep(NA, nrow(df)))
    out <- df[[cols[1]]]
    if (length(cols) > 1) {
      for (cc in cols[-1]) out <- dplyr::coalesce(out, df[[cc]])
    }
    out
  }

  collapse_stage <- function(x) {
    x <- toupper(as.character(x))
    dplyr::case_when(
      stringr::str_detect(x, "STAGE I[^V]|STAGE IA|STAGE IB|^I$|^IA$|^IB$") ~ "I",
      stringr::str_detect(x, "STAGE II|^II$|^IIA$|^IIB$") ~ "II",
      stringr::str_detect(x, "STAGE III|^III$|^IIIA$|^IIIB$|^IIIC$") ~ "III",
      stringr::str_detect(x, "STAGE IV|^IV$|^IVA$|^IVB$") ~ "IV",
      TRUE ~ NA_character_
    )
  }

  survdat <- scores2 |>
    dplyr::filter(.data$tissue == "tumor") |>
    dplyr::mutate(
      vital_status_tmp = as.character(coalesce_cols(dplyr::cur_data_all(), c("vital_status.y", "vital_status.x", "vital_status"))),
      days_to_death_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("days_to_death.y", "days_to_death.x", "days_to_death")))),
      days_to_follow_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("days_to_last_follow_up.y", "days_to_last_follow_up.x", "days_to_last_follow_up")))),
      age_tmp = suppressWarnings(as.numeric(coalesce_cols(dplyr::cur_data_all(), c("age_at_index.y", "age_at_index.x", "age_at_index")))),
      gender_tmp = as.factor(coalesce_cols(dplyr::cur_data_all(), c("gender.y", "gender.x", "gender", "sex"))),
      stage_raw = as.character(coalesce_cols(dplyr::cur_data_all(), c("ajcc_pathologic_stage.y", "ajcc_pathologic_stage.x", "ajcc_pathologic_stage"))),
      stage_collapsed = factor(collapse_stage(stage_raw), levels = c("I", "II", "III", "IV")),
      time_days = ifelse(!is.na(days_to_death_tmp), days_to_death_tmp, days_to_follow_tmp),
      event = ifelse(stringr::str_detect(vital_status_tmp, stringr::regex("dead", ignore_case = TRUE)), 1, 0),
      time_months = time_days / 30.44
    ) |>
    dplyr::filter(is.finite(.data$time_months), .data$time_months > 0)

  readr::write_tsv(survdat, "results/tables/tcga_lihc_survival_model_input.tsv")

  if (nrow(survdat) < 30 || length(unique(survdat$event)) < 2) {
    message("[06] Skipping Cox models: insufficient survival events or samples.")
    return(tibble::tibble())
  }

  has_age <- sum(!is.na(survdat$age_tmp)) >= 100
  has_gender <- nlevels(droplevels(survdat$gender_tmp)) > 1
  has_stage <- nlevels(droplevels(survdat$stage_collapsed)) > 1

  cov_age_sex <- character()
  if (has_age) cov_age_sex <- c(cov_age_sex, "age_tmp")
  if (has_gender) cov_age_sex <- c(cov_age_sex, "gender_tmp")

  cov_age_sex_stage <- cov_age_sex
  if (has_stage) cov_age_sex_stage <- c(cov_age_sex_stage, "stage_collapsed")

  fit_one <- function(score_name, covars = character(), model_name = "score_only") {
    rhs <- paste0("scale(", score_name, ")")
    if (length(covars) > 0) rhs <- paste(rhs, paste(covars, collapse = " + "), sep = " + ")
    f <- stats::as.formula(paste("survival::Surv(time_months, event) ~", rhs))
    fit <- tryCatch(
      survival::coxph(f, data = survdat),
      error = function(e) {
        message("[06] Cox model failed for ", score_name, " / ", model_name, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      dplyr::mutate(
        model = model_name,
        score = score_name,
        formula = paste(deparse(f), collapse = " "),
        n = fit$n,
        events = fit$nevent,
        .before = 1
      )
  }

  scores_to_test <- c("ProlifHubScore", "HepLossScore", "HCCStateScore")
  out <- dplyr::bind_rows(
    lapply(scores_to_test, fit_one, covars = character(), model_name = "score_only"),
    lapply(scores_to_test, fit_one, covars = cov_age_sex, model_name = "age_sex_adjusted"),
    lapply(scores_to_test, fit_one, covars = cov_age_sex_stage, model_name = "age_sex_stage_adjusted")
  )
  out
}

# -----------------------------
# Download and prepare TCGA-LIHC
# -----------------------------
query <- TCGAbiolinks::GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
TCGAbiolinks::GDCdownload(query)
se <- TCGAbiolinks::GDCprepare(query)

expr_counts <- SummarizedExperiment::assay(se)
expr_counts <- collapse_counts_to_symbols(expr_counts, SummarizedExperiment::rowData(se))
expr_log <- log2(edgeR::cpm(expr_counts, log = FALSE) + 1)
rownames(expr_log) <- toupper(rownames(expr_log))

meta <- make_coldata_tibble(se) |>
  dplyr::mutate(tissue = infer_tcga_tissue(dplyr::cur_data_all()))

# Keep only samples that can be interpreted as tumor or normal/non-tumor.
meta <- meta |>
  dplyr::filter(.data$sample_id %in% colnames(expr_log)) |>
  dplyr::arrange(match(.data$sample_id, colnames(expr_log)))
expr_log <- expr_log[, meta$sample_id, drop = FALSE]

label_audit <- meta |>
  dplyr::count(.data$tissue, name = "n") |>
  dplyr::mutate(dataset = "TCGA-LIHC", .before = 1)
readr::write_tsv(label_audit, "results/tables/tcga_lihc_tissue_label_audit.tsv")

if (!all(c("tumor", "non_tumor") %in% unique(meta$tissue))) {
  message("[06] Warning: TCGA-LIHC did not contain both tumor and non_tumor labels after parsing. Continuing with available labels.")
}

scores <- score_modules(expr_log, panels$ProlifHub, panels$HepLoss) |>
  dplyr::left_join(meta, by = "sample_id")
readr::write_tsv(scores, "results/tables/tcga_lihc_module_scores.tsv")

# Clinical metadata
clin <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LIHC", type = "clinical") |>
  tibble::as_tibble()
readr::write_tsv(clin, "results/tables/tcga_lihc_clinical_raw.tsv")
clin2 <- clean_clinical_for_merge(clin)

scores2 <- scores |>
  dplyr::mutate(patient = substr(.data$sample_id, 1, 12)) |>
  dplyr::left_join(clin2, by = "patient")
readr::write_tsv(scores2, "results/tables/tcga_lihc_module_scores_with_clinical.tsv")

assoc <- fit_lm_associations(scores2)
readr::write_tsv(assoc, "results/tables/tcga_lihc_clinicopathologic_associations.tsv")

cox <- fit_survival_models(scores2)
readr::write_tsv(cox, "results/tables/tcga_lihc_survival_cox_models.tsv")

message("[06] TCGA-LIHC pipeline completed. Outputs written to results/tables/.")
