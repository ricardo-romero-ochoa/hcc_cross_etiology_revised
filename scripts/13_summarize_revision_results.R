# ============================================================
# 13_summarize_revision_results.R
# Consolidates reviewer-facing summaries after the full pipeline.
# ============================================================

source("R/_shared.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

ensure_dirs()
dir.create("results/revision_summary", recursive = TRUE, showWarnings = FALSE)

# File inventory
expected <- c(
  "results/tables/tissue_label_audit.tsv",
  "results/tables/gene_symbol_overlap_audit.tsv",
  "results/tables/discovery_hallmark_combined.tsv",
  "results/tables/meta_analysis_with_heterogeneity.tsv",
  "results/tables/module_gene_panels_top20.tsv",
  "results/tables/geo_validation_summary.tsv",
  "results/tables/module_size_robustness.tsv",
  "results/tables/HBV_INJURY_topN_extended_summary.tsv",
  "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv",
  "results/tables/tcga_lihc_survival_cox_models.tsv"
)
readr::write_tsv(tibble(file = expected, exists = file.exists(expected)), "results/revision_summary/00_file_inventory.tsv")

# Meta-analysis summary
if (file.exists("results/tables/meta_analysis_with_heterogeneity.tsv")) {
  meta <- readr::read_tsv("results/tables/meta_analysis_with_heterogeneity.tsv", show_col_types = FALSE)
  meta_summary <- meta |>
    filter(meta_status == "ok") |>
    summarise(
      meta_analysis_rows = nrow(meta),
      estimable_genes = n(),
      non_estimable_genes = sum(meta$meta_status == "not_estimable", na.rm = TRUE),
      genes_FE_FDR_lt_0_05 = sum(meta_FDR_FE < 0.05, na.rm = TRUE),
      median_I2 = median(I2, na.rm = TRUE),
      mean_I2 = mean(I2, na.rm = TRUE),
      genes_I2_gt_50 = sum(I2 > 50, na.rm = TRUE),
      median_tau2_REML = median(tau2_REML, na.rm = TRUE),
      mean_tau2_REML = mean(tau2_REML, na.rm = TRUE),
      genes_tau2_REML_gt_0 = sum(tau2_REML > 0, na.rm = TRUE)
    )
  readr::write_tsv(meta_summary, "results/revision_summary/01_meta_analysis_summary.tsv")
}

# HBV injury coefficient summary
if (file.exists("results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv")) {
  x <- readr::read_tsv("results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv", show_col_types = FALSE)
  readr::write_tsv(x, "results/revision_summary/02_hbv_injury_adjusted_models.tsv")
}

# TCGA score terms only
if (file.exists("results/tables/tcga_lihc_survival_cox_models.tsv")) {
  cox <- readr::read_tsv("results/tables/tcga_lihc_survival_cox_models.tsv", show_col_types = FALSE)
  cox_terms <- cox |>
    filter(stringr::str_detect(term, "ProlifHubScore|HepLossScore|HCCStateScore"))
  readr::write_tsv(cox_terms, "results/revision_summary/03_tcga_survival_score_terms.tsv")
}

message("[13] Wrote revision summaries to results/revision_summary/.")
