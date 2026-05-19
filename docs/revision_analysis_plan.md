# Revision analysis plan mapped to reviewer concerns

## Reviewer concern: limited validation
Implemented in `scripts/05_validate_geo_cohorts.R` using GSE14520 plus additional GEO cohorts listed in `config/datasets.yml`. Outputs: `results/tables/geo_validation_summary.tsv` and `results/figures/Fig5_multicohort_validation.pdf`.

## Reviewer concern: unclear score construction
Implemented in `scripts/04_meta_modules.R` and `R/_shared.R`. The modules are generated from conserved genes with concordant direction, cohort-level FDR < 0.05, and meta-analysis FDR < 0.05. Scores are mean gene-wise z-scores. `HCCStateScore = ProlifHubScore - HepLossScore`.

## Reviewer concern: overfitting / arbitrary top 20 genes
Implemented in `scripts/10_module_size_sensitivity.R`. The pipeline rebuilds modules for top 10, 15, 20, 30, and 50 genes and evaluates score performance across cohorts.

## Reviewer concern: missing meta-analysis heterogeneity
Implemented in `scripts/04_meta_modules.R`. Outputs include I2, tau2, Q, and Q-test p-value, plus fixed-effect and REML random-effects sensitivity estimates.

## Reviewer concern: tissue composition / immune infiltration
Implemented in two complementary ways:

1. `scripts/07_estimate_adjustment.R` runs ESTIMATE when possible.
2. `scripts/08_cibersortx_export_GSE121248.R` and `scripts/09_cibersortx_adjusted_regression_GSE121248.R` export GSE121248 for CIBERSORTx and fit adjusted regression models after CIBERSORTx results are returned.

## Reviewer concern: HBV_INJURY gene list absent
Implemented in `scripts/03b_derive_hepatitis_axes.R`, which writes `results/tables/HBV_INJURY_gene_set.tsv` with the gene list and derivation statistics.

## Reviewer concern: TCGA underused
Implemented in `scripts/06_tcga_lihc_pipeline.R`, which computes module scores, clinicopathologic associations, and Cox survival models in TCGA-LIHC.

## Reviewer concern: figures and panel labels
Implemented in `scripts/11_make_revision_figures.R`. These figures are meant as clean templates; additional manuscript styling can be applied after confirming all results.
