source("R/_shared.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
  library(broom)
})

ensure_dirs()
cfg <- read_config()
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

find_hallmark_col <- function(df, pattern) {
  hits <- grep(pattern, colnames(df), value = TRUE, ignore.case = TRUE)

  if (length(hits) == 0) {
    stop(
      "Could not find Hallmark column matching pattern: ", pattern,
      "\nAvailable columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }

  hits[1]
}

load_hallmark_scores_wide <- function(path) {
  if (!file.exists(path)) {
    stop("Missing Hallmark score file: ", path,
         "\nRun scripts/03_discovery_hallmark_gsva.R first.")
  }

  obj <- read_rds(path)

  # Accept either:
  # 1. list(scores = matrix/pathway x samples)
  # 2. matrix directly
  # 3. data.frame/tibble
  hall <- if (is.list(obj) && "scores" %in% names(obj)) obj$scores else obj

  if (is.matrix(hall)) {
    # Expected: rows = pathways, columns = samples
    hall_df <- as.data.frame(t(hall)) |>
      rownames_to_column("sample_id") |>
      as_tibble()
  } else {
    hall_df <- as_tibble(hall)

    # If already sample-wise, keep as is.
    if (!("sample_id" %in% colnames(hall_df))) {
      first_col <- colnames(hall_df)[1]

      # Assume rows are pathways and columns are samples.
      hall_df <- hall_df |>
        rename(pathway = all_of(first_col)) |>
        pivot_longer(
          cols = -pathway,
          names_to = "sample_id",
          values_to = "score"
        ) |>
        pivot_wider(
          names_from = pathway,
          values_from = score
        )
    }
  }

  e2f_col <- find_hallmark_col(hall_df, "E2F.*TARGETS")
  g2m_col <- find_hallmark_col(hall_df, "G2M.*CHECKPOINT")

  hall_df |>
    transmute(
      sample_id = as.character(sample_id),
      E2F = as.numeric(.data[[e2f_col]]),
      G2M = as.numeric(.data[[g2m_col]])
    )
}

safe_lm <- function(formula, data, model_name) {
  dat <- model.frame(formula, data = data, na.action = na.omit)

  if (nrow(dat) < 10) {
    warning("Skipping model ", model_name, ": fewer than 10 complete samples.")
    return(NULL)
  }

  fit <- tryCatch(
    lm(formula, data = data),
    error = function(e) {
      warning("Skipping model ", model_name, ": ", conditionMessage(e))
      NULL
    }
  )

  fit
}

read_injury_genes <- function(path) {
  if (!file.exists(path)) {
    stop("Missing injury gene-set file: ", path)
  }

  tbl <- readr::read_tsv(path, show_col_types = FALSE)

  if (!("gene" %in% colnames(tbl))) {
    stop("Injury gene-set file must contain a 'gene' column: ", path)
  }

  genes <- tbl$gene |>
    as.character() |>
    unique() |>
    na.omit()

  genes[genes != ""]
}

# ------------------------------------------------------------
# Optional: create top-N/extended HBV_INJURY sets if missing
# ------------------------------------------------------------
# The script will use existing files if present. If they are absent and
# results/tables/GSE83148_HBV_INJURY_derivation_full.tsv exists, it will
# generate top-N gene-set files automatically.
# ------------------------------------------------------------

create_hbv_injury_topn_sets_if_needed <- function(
  deriv_path = "results/tables/GSE83148_HBV_INJURY_derivation_full.tsv",
  top_ns = c(200, 500, 1000, 2000, 5000)
) {
  expected <- c(
    paste0("results/tables/HBV_INJURY_TOP_", top_ns, "_gene_set.tsv"),
    "results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv"
  )

  if (all(file.exists(expected))) {
    return(invisible(TRUE))
  }

  if (!file.exists(deriv_path)) {
    warning(
      "Some HBV_INJURY top-N files are missing and derivation file is unavailable: ",
      deriv_path,
      "\nRun scripts/03c_make_hbv_injury_gene_set.R and the top-N gene-set generation step first."
    )
    return(invisible(FALSE))
  }

  message("[09] Creating missing HBV_INJURY top-N/extended gene-set files from ", deriv_path)

  tt <- readr::read_tsv(deriv_path, show_col_types = FALSE)

  required <- c("gene", "beta_injury_index", "p_value", "FDR")
  missing <- setdiff(required, colnames(tt))
  if (length(missing) > 0) {
    stop("Derivation table is missing columns: ", paste(missing, collapse = ", "))
  }

  if (!("t" %in% colnames(tt))) {
    tt$t <- NA_real_
  }

  ranked <- tt |>
    filter(beta_injury_index > 0, FDR < 0.10) |>
    arrange(FDR, desc(abs(t)), desc(beta_injury_index)) |>
    mutate(
      rank = row_number(),
      selection_rule = "positive_beta_FDR_lt_0.10_ranked_by_FDR_and_moderated_t"
    )

  if (nrow(ranked) < 10) {
    stop("Too few positive FDR<0.10 injury-associated genes in derivation table: ", nrow(ranked))
  }

  extended <- ranked |>
    mutate(injury_set = paste0("HBV_INJURY_EXTENDED_", nrow(ranked))) |>
    select(injury_set, rank, gene, beta_injury_index, t, p_value, FDR, selection_rule)

  extended_path <- paste0("results/tables/HBV_INJURY_EXTENDED_", nrow(ranked), "_gene_set.tsv")
  readr::write_tsv(extended, extended_path)

  # Keep the expected 7792 filename for compatibility if this is the observed extended size.
  # If the size differs, still create the conventional path containing the true extended set.
  readr::write_tsv(extended, "results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv")

  for (n in top_ns) {
    nn <- min(n, nrow(ranked))
    out <- ranked |>
      slice_head(n = nn) |>
      mutate(
        injury_set = paste0("HBV_INJURY_TOP_", n),
        selection_rule = paste0("top_", n, "_positive_beta_FDR_lt_0.10_ranked")
      ) |>
      select(injury_set, rank, gene, beta_injury_index, t, p_value, FDR, selection_rule)

    readr::write_tsv(out, paste0("results/tables/HBV_INJURY_TOP_", n, "_gene_set.tsv"))
  }

  summary_tbl <- bind_rows(
    purrr::map_dfr(top_ns, function(n) {
      path <- paste0("results/tables/HBV_INJURY_TOP_", n, "_gene_set.tsv")
      readr::read_tsv(path, show_col_types = FALSE)
    }),
    extended
  ) |>
    group_by(injury_set, selection_rule) |>
    summarise(
      n_genes = n(),
      min_FDR = min(FDR, na.rm = TRUE),
      max_FDR = max(FDR, na.rm = TRUE),
      median_beta = median(beta_injury_index, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_tsv(summary_tbl, "results/tables/HBV_INJURY_topN_extended_summary.tsv")
  invisible(TRUE)
}

create_hbv_injury_topn_sets_if_needed()

# ------------------------------------------------------------
# Load curated GSE121248 expression and metadata
# ------------------------------------------------------------

obj <- read_rds("data/processed/GSE121248_curated.rds")

expr <- obj$expr

meta <- obj$metadata |>
  mutate(sample_id = as.character(sample_id)) |>
  select(sample_id, tissue, everything())

if (!all(colnames(expr) %in% meta$sample_id)) {
  warning("Some expression samples are missing from metadata.")
}

# ------------------------------------------------------------
# Load Hallmark E2F/G2M scores robustly
# ------------------------------------------------------------

hall_df <- load_hallmark_scores_wide(
  "data/processed/GSE121248_hallmark_scores.rds"
)

# ------------------------------------------------------------
# Import CIBERSORTx result
# ------------------------------------------------------------

if (is.null(cfg$cibersortx$results_in) || !file.exists(cfg$cibersortx$results_in)) {
  stop(
    "CIBERSORTx results not found: ", cfg$cibersortx$results_in,
    "\nRun CIBERSORTx externally first and save the result to the configured path."
  )
}

cs <- readr::read_tsv(cfg$cibersortx$results_in, show_col_types = FALSE)

# CIBERSORTx usually uses "Mixture"; some exported/edited files may already use "sample_id".
sample_col <- dplyr::case_when(
  "sample_id" %in% colnames(cs) ~ "sample_id",
  "Mixture" %in% colnames(cs) ~ "Mixture",
  "Mixture ID" %in% colnames(cs) ~ "Mixture ID",
  TRUE ~ colnames(cs)[1]
)

if (!sample_col %in% colnames(cs)) {
  stop(
    "Could not identify the CIBERSORTx sample column. Columns found:\n",
    paste(colnames(cs), collapse = ", ")
  )
}

if (sample_col != "sample_id") {
  colnames(cs)[colnames(cs) == sample_col] <- "sample_id"
}

cs <- cs |>
  mutate(sample_id = as.character(sample_id))

drop_cols <- c(
  "P-value", "P.value", "P_value",
  "Correlation",
  "RMSE",
  "Absolute score",
  "Absolute_score"
)

frac_cols <- setdiff(colnames(cs), c("sample_id", drop_cols))

if (length(frac_cols) < 2) {
  stop(
    "Too few CIBERSORTx fraction columns detected.\n",
    "Columns found:\n",
    paste(colnames(cs), collapse = ", ")
  )
}

frac <- cs |>
  select(sample_id, all_of(frac_cols)) |>
  mutate(across(-sample_id, ~ suppressWarnings(as.numeric(.x))))

# Keep only expression samples in the CIBERSORTx fraction table.
frac <- frac |>
  filter(sample_id %in% colnames(expr))

if (nrow(frac) < 10) {
  stop("Too few CIBERSORTx samples overlap with GSE121248 expression matrix.")
}

# ------------------------------------------------------------
# PC model to address compositional collinearity
# ------------------------------------------------------------

frac_mat <- frac |>
  column_to_rownames("sample_id") |>
  as.matrix()

# Remove all-NA or zero-variance fractions.
keep_frac <- apply(frac_mat, 2, function(x) {
  sum(is.finite(x)) >= 5 && stats::var(x, na.rm = TRUE) > 0
})

frac_mat <- frac_mat[, keep_frac, drop = FALSE]

if (ncol(frac_mat) < 2) {
  stop("Too few usable CIBERSORTx fractions after removing zero-variance columns.")
}

# Replace remaining NA values with column medians.
for (j in seq_len(ncol(frac_mat))) {
  x <- frac_mat[, j]
  if (anyNA(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    frac_mat[, j] <- x
  }
}

frac_clr <- clr_transform(frac_mat)

n_pcs <- min(cfg$cibersortx$n_pcs %||% 2, ncol(frac_clr), nrow(frac_clr) - 1)

pc <- prcomp(frac_clr, center = TRUE, scale. = TRUE)

pc_df <- as_tibble(
  pc$x[, seq_len(n_pcs), drop = FALSE],
  rownames = "sample_id"
) |>
  rename_with(
    ~ paste0("CIBERSORTx_PC", seq_along(.x)),
    -sample_id
  )

frac_df <- as_tibble(frac_mat, rownames = "sample_id")

# Fraction summary does not depend on injury gene set.
fraction_summary <- frac_df |>
  pivot_longer(
    cols = -sample_id,
    names_to = "cell_fraction",
    values_to = "fraction"
  ) |>
  left_join(meta |> select(sample_id, tissue), by = "sample_id") |>
  group_by(tissue, cell_fraction) |>
  summarise(
    n = sum(!is.na(fraction)),
    mean = mean(fraction, na.rm = TRUE),
    median = median(fraction, na.rm = TRUE),
    sd = sd(fraction, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  fraction_summary,
  "results/tables/gse121248_cibersortx_cell_fraction_summary.tsv"
)

# ------------------------------------------------------------
# Run one HBV_INJURY gene-set regression
# ------------------------------------------------------------

run_hbv_injury_regression <- function(
  injury_path,
  injury_set_label,
  analysis_role,
  expr,
  meta,
  hall_df,
  pc_df,
  frac_df,
  frac_mat
) {
  injury_genes <- read_injury_genes(injury_path)
  overlap_genes <- intersect(injury_genes, rownames(expr))
  overlap_n <- length(overlap_genes)

  if (overlap_n < 5) {
    stop(
      "Too few genes overlap expression matrix for ", injury_set_label,
      ". Overlap = ", overlap_n
    )
  }

  hbv_injury <- tibble(
    sample_id = colnames(expr),
    HBV_INJURY = score_module(expr, injury_genes, min_genes = 5)
  )

  base_dat <- meta |>
    select(sample_id, tissue) |>
    left_join(hbv_injury, by = "sample_id") |>
    left_join(hall_df, by = "sample_id") |>
    filter(tissue %in% c("tumor", "non_tumor")) |>
    mutate(
      tissue = factor(tissue, levels = c("non_tumor", "tumor"))
    ) |>
    filter(
      !is.na(HBV_INJURY),
      !is.na(E2F),
      !is.na(G2M)
    )

  if (nrow(base_dat) < 10) {
    stop("Too few complete samples for ", injury_set_label, " after merging metadata/HBV_INJURY/E2F/G2M.")
  }

  dat <- base_dat |>
    left_join(pc_df, by = "sample_id") |>
    left_join(frac_df, by = "sample_id")

  pc_terms <- colnames(pc_df)[colnames(pc_df) != "sample_id"]

  models <- list()

  models$unadjusted <- safe_lm(
    HBV_INJURY ~ tissue,
    data = dat,
    model_name = paste(injury_set_label, "unadjusted")
  )

  models$proliferation_adjusted <- safe_lm(
    HBV_INJURY ~ tissue + E2F + G2M,
    data = dat,
    model_name = paste(injury_set_label, "proliferation_adjusted")
  )

  if (length(pc_terms) > 0) {
    pc_formula <- as.formula(
      paste(
        "HBV_INJURY ~ tissue + E2F + G2M +",
        paste(pc_terms, collapse = " + ")
      )
    )

    models$proliferation_cibersortx_pc_adjusted <- safe_lm(
      pc_formula,
      data = dat,
      model_name = paste(injury_set_label, "proliferation_cibersortx_pc_adjusted")
    )
  }

  # Selected-fraction model: keep top variable fractions and use syntactic names.
  frac_vars <- sort(
    apply(frac_mat, 2, var, na.rm = TRUE),
    decreasing = TRUE
  )

  selected_raw <- names(frac_vars)[seq_len(min(8, length(frac_vars)))]

  dat_fraction <- dat
  colnames(dat_fraction) <- make.names(colnames(dat_fraction), unique = TRUE)

  selected <- make.names(selected_raw, unique = TRUE)
  selected <- selected[selected %in% colnames(dat_fraction)]

  if (length(selected) > 0) {
    selected_formula <- as.formula(
      paste(
        "HBV_INJURY ~ tissue + E2F + G2M +",
        paste(selected, collapse = " + ")
      )
    )

    models$proliferation_selected_fraction_adjusted <- safe_lm(
      selected_formula,
      data = dat_fraction,
      model_name = paste(injury_set_label, "proliferation_selected_fraction_adjusted")
    )
  }

  models <- models[!vapply(models, is.null, logical(1))]

  if (length(models) == 0) {
    stop("No regression models could be fitted for ", injury_set_label)
  }

  coef_tab <- imap_dfr(
    models,
    function(fit, nm) {
      broom::tidy(fit, conf.int = TRUE) |>
        mutate(
          analysis_role = analysis_role,
          injury_set = injury_set_label,
          injury_gene_file = injury_path,
          n_input_genes = length(injury_genes),
          n_overlap_genes = overlap_n,
          model = nm,
          .before = 1
        )
    }
  )

  tissue_only <- coef_tab |>
    filter(term == "tissuetumor")

  unadj <- tissue_only$estimate[tissue_only$model == "unadjusted"][1]

  tissue_only <- tissue_only |>
    mutate(
      percent_retained_vs_unadjusted = 100 * estimate / unadj
    )

  input_dat <- dat |>
    mutate(
      analysis_role = analysis_role,
      injury_set = injury_set_label,
      injury_gene_file = injury_path,
      n_input_genes = length(injury_genes),
      n_overlap_genes = overlap_n,
      .before = 1
    )

  list(
    full = coef_tab,
    tissue_only = tissue_only,
    input = input_dat
  )
}

# ------------------------------------------------------------
# Gene sets to evaluate
# ------------------------------------------------------------
# Interpretation:
# - TOP_200/500/1000/2000/5000 are compact top-N sensitivity sets.
# - EXTENDED_7792 is the full positive FDR<0.10 injury program.
# Use the smallest top-N set that is stable across adjustment models as a
# compact primary candidate; report EXTENDED_7792 as extended sensitivity.
# ------------------------------------------------------------

injury_sets <- tibble::tribble(
  ~injury_set_label,          ~injury_path,                                           ~analysis_role,
  "HBV_INJURY_TOP_200",       "results/tables/HBV_INJURY_TOP_200_gene_set.tsv",       "topN_sensitivity",
  "HBV_INJURY_TOP_500",       "results/tables/HBV_INJURY_TOP_500_gene_set.tsv",       "topN_sensitivity",
  "HBV_INJURY_TOP_1000",      "results/tables/HBV_INJURY_TOP_1000_gene_set.tsv",      "topN_sensitivity",
  "HBV_INJURY_TOP_2000",      "results/tables/HBV_INJURY_TOP_2000_gene_set.tsv",      "compact_primary_candidate",
  "HBV_INJURY_TOP_5000",      "results/tables/HBV_INJURY_TOP_5000_gene_set.tsv",      "topN_sensitivity",
  "HBV_INJURY_EXTENDED_7792", "results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv", "extended_sensitivity"
) |>
  filter(file.exists(injury_path))

if (nrow(injury_sets) == 0) {
  stop("No HBV_INJURY top-N/extended gene-set files found in results/tables/.")
}

message("[09] Running HBV_INJURY regression for ", nrow(injury_sets), " gene-set definitions:")
print(injury_sets)

res <- purrr::pmap(
  injury_sets,
  function(injury_set_label, injury_path, analysis_role) {
    run_hbv_injury_regression(
      injury_path = injury_path,
      injury_set_label = injury_set_label,
      analysis_role = analysis_role,
      expr = expr,
      meta = meta,
      hall_df = hall_df,
      pc_df = pc_df,
      frac_df = frac_df,
      frac_mat = frac_mat
    )
  }
)

coef_full <- purrr::map_dfr(res, "full")
coef_tissue <- purrr::map_dfr(res, "tissue_only")
reg_inputs <- purrr::map_dfr(res, "input")

readr::write_tsv(
  coef_full,
  "results/tables/gse121248_hbv_injury_topN_extended_regression_full.tsv"
)

readr::write_tsv(
  coef_tissue,
  "results/tables/gse121248_hbv_injury_topN_extended_regression_tissue.tsv"
)

readr::write_tsv(
  reg_inputs,
  "results/tables/gse121248_hbv_injury_topN_extended_regression_input.tsv"
)

manuscript_summary <- coef_tissue |>
  transmute(
    analysis_role,
    injury_set,
    n_input_genes,
    n_overlap_genes,
    model,
    tumor_coefficient = estimate,
    percent_retained_vs_unadjusted,
    p_value = p.value,
    ci_low = conf.low,
    ci_high = conf.high
  )

readr::write_tsv(
  manuscript_summary,
  "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv"
)

# ------------------------------------------------------------
# Backward-compatible outputs
# ------------------------------------------------------------
# For older figure/table scripts, write conventional CIBERSORTx output files.
# Preferred source:
#   1. HBV_INJURY_TOP_2000 if available; otherwise
#   2. HBV_INJURY_EXTENDED_7792; otherwise
#   3. first evaluated set.
# ------------------------------------------------------------

preferred_set <- dplyr::case_when(
  "HBV_INJURY_TOP_2000" %in% manuscript_summary$injury_set ~ "HBV_INJURY_TOP_2000",
  "HBV_INJURY_EXTENDED_7792" %in% manuscript_summary$injury_set ~ "HBV_INJURY_EXTENDED_7792",
  TRUE ~ manuscript_summary$injury_set[1]
)

coef_full_pref <- coef_full |>
  filter(injury_set == preferred_set) |>
  select(-analysis_role, -injury_set, -injury_gene_file, -n_input_genes, -n_overlap_genes)

coef_tissue_pref <- coef_tissue |>
  filter(injury_set == preferred_set) |>
  select(model, estimate, std.error, statistic, p.value, conf.low, conf.high)

reg_input_pref <- reg_inputs |>
  filter(injury_set == preferred_set) |>
  select(-analysis_role, -injury_set, -injury_gene_file, -n_input_genes, -n_overlap_genes)

readr::write_tsv(
  coef_full_pref,
  "results/tables/gse121248_cibersortx_adjusted_regression_full.tsv"
)

readr::write_tsv(
  coef_tissue_pref,
  "results/tables/gse121248_cibersortx_adjusted_regression.tsv"
)

readr::write_tsv(
  reg_input_pref,
  "results/tables/gse121248_cibersortx_regression_input.tsv"
)

message("[09] Preferred backward-compatible set: ", preferred_set)
message("[09] Wrote:")
message("  results/tables/gse121248_hbv_injury_topN_extended_regression_full.tsv")
message("  results/tables/gse121248_hbv_injury_topN_extended_regression_tissue.tsv")
message("  results/tables/gse121248_hbv_injury_topN_extended_regression_input.tsv")
message("  results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv")
message("  results/tables/gse121248_cibersortx_adjusted_regression_full.tsv")
message("  results/tables/gse121248_cibersortx_adjusted_regression.tsv")
message("  results/tables/gse121248_cibersortx_regression_input.tsv")
message("  results/tables/gse121248_cibersortx_cell_fraction_summary.tsv")
