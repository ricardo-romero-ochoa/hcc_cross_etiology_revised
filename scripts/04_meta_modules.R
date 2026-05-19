source("R/_shared.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(metafor)
})

ensure_dirs()
cfg <- read_config()

message("[04] Running gene-level limma contrasts for discovery cohorts")
res <- list()
for (acc in c("GSE121248", "GSE41804")) {
  obj_path <- paste0("data/processed/", acc, "_curated.rds")
  if (!file.exists(obj_path)) stop("Missing curated object: ", obj_path, ". Run script 02 first.")
  obj <- read_rds(obj_path)
  de <- limma_tumor_contrast(obj$expr, obj$metadata)$table |>
    mutate(dataset = acc)
  readr::write_tsv(de, paste0("results/tables/", acc, "_gene_limma.tsv"))
  res[[acc]] <- de
}

required_cols <- c("gene", "logFC", "t", "p_value", "FDR")
for (acc in names(res)) {
  missing_cols <- setdiff(required_cols, names(res[[acc]]))
  if (length(missing_cols) > 0) {
    stop("Limma table for ", acc, " is missing required columns: ", paste(missing_cols, collapse = ", "),
         ". Columns found: ", paste(names(res[[acc]]), collapse = ", "))
  }
}

wide <- full_join(
  res[["GSE121248"]] |>
    select(gene, logFC_121248 = logFC, t_121248 = t, p_121248 = p_value, FDR_121248 = FDR),
  res[["GSE41804"]] |>
    select(gene, logFC_41804 = logFC, t_41804 = t, p_41804 = p_value, FDR_41804 = FDR),
  by = "gene"
)

readr::write_tsv(wide, "results/tables/discovery_gene_limma_wide.tsv")

# Return a one-row tibble even when a gene cannot be meta-analyzed. This avoids
# dplyr errors such as `meta_p_FE not found` when many rows are invalid or no rows
# pass the SE checks.
empty_meta_row <- function(g) {
  tibble(
    gene = g,
    meta_logFC_FE = NA_real_,
    meta_se_FE = NA_real_,
    meta_z_FE = NA_real_,
    meta_p_FE = NA_real_,
    Q = NA_real_,
    Q_p = NA_real_,
    I2 = NA_real_,
    tau2_FE = NA_real_,
    meta_logFC_REML = NA_real_,
    meta_p_REML = NA_real_,
    tau2_REML = NA_real_,
    meta_status = "not_estimable"
  )
}

# Approximate standard errors from limma moderated t-statistics: SE = logFC / t.
# This is acceptable for ranking and heterogeneity reporting but should be described
# as an approximation in the methods.
meta_one <- function(g, l1, t1, l2, t2) {
  if (any(is.na(c(l1, t1, l2, t2))) || any(c(t1, t2) == 0)) {
    return(empty_meta_row(g))
  }
  sei <- abs(c(l1 / t1, l2 / t2))
  yi <- c(l1, l2)
  if (any(!is.finite(sei)) || any(sei <= 0) || any(!is.finite(yi))) {
    return(empty_meta_row(g))
  }
  fe <- tryCatch(metafor::rma.uni(yi = yi, sei = sei, method = "FE"), error = function(e) e)
  re <- tryCatch(metafor::rma.uni(yi = yi, sei = sei, method = "REML"), error = function(e) e)
  if (inherits(fe, "error")) return(empty_meta_row(g))

  tibble(
    gene = g,
    meta_logFC_FE = as.numeric(fe$b),
    meta_se_FE = as.numeric(fe$se),
    meta_z_FE = as.numeric(fe$zval),
    meta_p_FE = as.numeric(fe$pval),
    Q = as.numeric(fe$QE),
    Q_p = as.numeric(fe$QEp),
    I2 = as.numeric(fe$I2),
    tau2_FE = as.numeric(fe$tau2),
    meta_logFC_REML = if (!inherits(re, "error")) as.numeric(re$b) else NA_real_,
    meta_p_REML = if (!inherits(re, "error")) as.numeric(re$pval) else NA_real_,
    tau2_REML = if (!inherits(re, "error")) as.numeric(re$tau2) else NA_real_,
    meta_status = "ok"
  )
}

message("[04] Running fixed-effect meta-analysis and REML sensitivity analysis")
meta_raw <- pmap_dfr(
  list(wide$gene, wide$logFC_121248, wide$t_121248, wide$logFC_41804, wide$t_41804),
  meta_one
)

if (!"meta_p_FE" %in% names(meta_raw)) {
  stop("Internal error: meta-analysis output lacks meta_p_FE. Columns found: ", paste(names(meta_raw), collapse = ", "))
}

meta <- meta_raw |>
  mutate(meta_FDR_FE = p.adjust(meta_p_FE, method = "BH")) |>
  left_join(wide, by = "gene") |>
  mutate(
    concordant_direction = is.finite(logFC_121248) & is.finite(logFC_41804) &
      sign(logFC_121248) == sign(logFC_41804),
    conserved = meta_status == "ok" & concordant_direction &
      FDR_121248 < cfg$fdr_cutoff & FDR_41804 < cfg$fdr_cutoff & meta_FDR_FE < cfg$fdr_cutoff,
    rank_sum_abs_t = rank(-abs(t_121248), ties.method = "average", na.last = "keep") +
      rank(-abs(t_41804), ties.method = "average", na.last = "keep")
  )

readr::write_tsv(meta, "results/tables/meta_analysis_with_heterogeneity.tsv")

n_ok <- sum(meta$meta_status == "ok", na.rm = TRUE)
n_cons <- sum(meta$conserved, na.rm = TRUE)
message("[04] Meta-analyzable genes: ", n_ok, "; conserved genes: ", n_cons)

if (n_cons == 0) {
  stop(
    "No conserved genes passed the configured filters. Inspect results/tables/discovery_gene_limma_wide.tsv ",
    "and results/tables/meta_analysis_with_heterogeneity.tsv. This usually means tissue labels are wrong, ",
    "gene mapping failed, or cfg$fdr_cutoff is too stringent for the current data."
  )
}

for (n in cfg$module_sizes) {
  up <- meta |>
    filter(conserved, meta_logFC_FE > 0) |>
    arrange(rank_sum_abs_t) |>
    slice_head(n = n) |>
    transmute(module = "ProlifHub", gene, rank = row_number(), module_size = n)

  down <- meta |>
    filter(conserved, meta_logFC_FE < 0) |>
    arrange(rank_sum_abs_t) |>
    slice_head(n = n) |>
    transmute(module = "HepLoss", gene, rank = row_number(), module_size = n)

  if (nrow(up) < n || nrow(down) < n) {
    warning("Requested top ", n, " genes, but available conserved genes were ProlifHub=", nrow(up),
            ", HepLoss=", nrow(down), ". Writing available genes only.")
  }

  readr::write_tsv(bind_rows(up, down), paste0("results/tables/module_gene_panels_top", n, ".tsv"))
}

primary_path <- paste0("results/tables/module_gene_panels_top", cfg$primary_module_size, ".tsv")
if (!file.exists(primary_path)) stop("Primary module file was not created: ", primary_path)
primary_dest <- "results/tables/module_gene_panels_top20.tsv"
if (normalizePath(primary_path, mustWork = FALSE) != normalizePath(primary_dest, mustWork = FALSE)) {
  file.copy(primary_path, primary_dest, overwrite = TRUE)
} else {
  message("[04] Primary module panel already written to ", primary_dest)
}

message("[04] Done. Outputs written to results/tables/meta_analysis_with_heterogeneity.tsv and module_gene_panels_top*.tsv")
