# ============================================================
# 12_make_manuscript_tables.R
# Writes compact manuscript and supplementary tables.
# ============================================================

source("R/_shared.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

ensure_dirs()
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)

# Table 1: gene-level meta-analysis summary and heterogeneity
if (file.exists("results/tables/meta_analysis_with_heterogeneity.tsv")) {
  meta <- readr::read_tsv("results/tables/meta_analysis_with_heterogeneity.tsv", show_col_types = FALSE)
  tbl1 <- meta |>
    filter(meta_status == "ok") |>
    summarise(
      `Meta-analysis rows` = nrow(meta),
      `Estimable genes` = n(),
      `Non-estimable genes` = sum(meta$meta_status == "not_estimable", na.rm = TRUE),
      `Genes with fixed-effect meta-analysis FDR < 0.05` = sum(meta_FDR_FE < 0.05, na.rm = TRUE),
      `Median I²` = median(I2, na.rm = TRUE),
      `Mean I²` = mean(I2, na.rm = TRUE),
      `Genes with I² > 50` = sum(I2 > 50, na.rm = TRUE),
      `Median τ² (REML)` = median(tau2_REML, na.rm = TRUE),
      `Mean τ² (REML)` = mean(tau2_REML, na.rm = TRUE),
      `Genes with τ² (REML) > 0` = sum(tau2_REML > 0, na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(everything(), names_to = "Metric", values_to = "Value")
  readr::write_csv(tbl1, "manuscript/tables/Table1_meta_analysis_summary.csv")

  conserved <- meta |>
    filter(conserved) |>
    arrange(rank_sum_abs_t) |>
    select(gene, meta_logFC_FE, meta_FDR_FE, I2, tau2_REML, Q_p, logFC_121248, FDR_121248, logFC_41804, FDR_41804)
  readr::write_csv(conserved, "manuscript/tables/SuppTable_conserved_meta_genes.csv")
}

# Table 2: HBV injury adjusted models
if (file.exists("results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv")) {
  hbv <- readr::read_tsv("results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv", show_col_types = FALSE) |>
    filter(injury_set %in% c("HBV_INJURY_TOP_2000", "HBV_INJURY_EXTENDED_7792")) |>
    mutate(
      `HBV injury set` = dplyr::recode(injury_set, HBV_INJURY_TOP_2000 = "TOP 2000", HBV_INJURY_EXTENDED_7792 = "EXTENDED 7792"),
      Role = dplyr::case_when(analysis_role == "compact_primary" ~ "Primary compact", TRUE ~ "Sensitivity"),
      Model = dplyr::recode(model,
        unadjusted = "Unadjusted",
        proliferation_adjusted = "E2F/G2M adjusted",
        proliferation_cibersortx_pc_adjusted = "E2F/G2M + CIBERSORTx PC",
        proliferation_selected_fraction_adjusted = "E2F/G2M + selected fractions"
      ),
      `Tumor coefficient` = round(tumor_coefficient, 3),
      `Effect retained vs unadjusted` = paste0(round(percent_retained_vs_unadjusted, 1), "%"),
      `p-value` = signif(p_value, 3),
      `95% CI` = paste0(round(ci_low, 4), "–", round(ci_high, 4))
    ) |>
    select(`HBV injury set`, Role, Model, `Tumor coefficient`, `Effect retained vs unadjusted`, `p-value`, `95% CI`)
  readr::write_csv(hbv, "manuscript/tables/Table2_HBV_injury_adjusted_models.csv")
}

# Table 3: TCGA-LIHC Cox models
if (file.exists("results/tables/tcga_lihc_survival_cox_models.tsv")) {
  cox <- readr::read_tsv("results/tables/tcga_lihc_survival_cox_models.tsv", show_col_types = FALSE) |>
    filter(stringr::str_detect(term, "ProlifHubScore|HepLossScore|HCCStateScore")) |>
    mutate(
      Model = dplyr::recode(model, score_only = "Score only", age_sex_adjusted = "Age/sex adjusted", age_sex_stage_adjusted = "Age/sex/stage adjusted"),
      Score = score,
      `HR per SD` = round(estimate, 3),
      `95% CI` = paste0(round(conf.low, 3), "–", round(conf.high, 3)),
      `p-value` = signif(p.value, 3)
    ) |>
    select(Model, Score, `HR per SD`, `95% CI`, `p-value`, n, events)
  readr::write_csv(cox, "manuscript/tables/Table3_TCGA_LIHC_Cox_models.csv")
}

# Validation and gene-set supplementary tables
if (file.exists("results/tables/geo_validation_summary.tsv")) {
  val <- readr::read_tsv("results/tables/geo_validation_summary.tsv", show_col_types = FALSE) |>
    mutate(across(where(is.numeric), ~round(.x, 4)))
  readr::write_csv(val, "manuscript/tables/SuppTable_multicohort_validation.csv")
}

if (file.exists("results/tables/HBV_INJURY_TOP_2000_gene_set.tsv")) {
  file.copy("results/tables/HBV_INJURY_TOP_2000_gene_set.tsv", "manuscript/tables/SuppTable_HBV_INJURY_TOP_2000_gene_set.tsv", overwrite = TRUE)
}
if (file.exists("results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv")) {
  file.copy("results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv", "manuscript/tables/SuppTable_HBV_INJURY_EXTENDED_7792_gene_set.tsv", overwrite = TRUE)
}

message("[12] Manuscript tables written to manuscript/tables/.")
