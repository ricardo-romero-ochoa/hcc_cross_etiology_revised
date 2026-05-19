[![DOI](https://zenodo.org/badge/1243934990.svg)](https://doi.org/10.5281/zenodo.20298126)


# HCC cross-etiology transcriptomic revision pipeline

This repository contains the complete reproducible R pipeline for the revised manuscript:

**Cross-etiology transcriptomic conservation in hepatocellular carcinoma reveals opposing proliferation and hepatocyte-loss programs validated across cohorts**

The workflow implements a cpmrehensive revision package: expanded GEO validation, explicit cohort curation, Hallmark GSVA, gene-level meta-analysis with heterogeneity, robust module-score construction, module-size sensitivity analysis, HBV injury-axis derivation, CIBERSORTx-adjusted regression, and TCGA-LIHC Cox modeling.

## Main analyses

The pipeline performs the following steps:

1. Curates GEO and TCGA-LIHC datasets.
2. Parses and audits tumor/non-tumor labels.
3. Maps microarray probes to HGNC symbols.
4. Scores MSigDB Hallmark pathways with GSVA/ssGSEA-compatible workflows.
5. Performs limma differential modeling in HBV-HCC and HCV-HCC discovery cohorts.
6. Performs inverse-variance fixed-effect meta-analysis and REML random-effects sensitivity analysis.
7. Reports heterogeneity using I² and REML-estimated τ².
8. Constructs ProlifHubScore, HepLossScore, and HCCStateScore.
9. Validates HCCStateScore across independent GEO HCC cohorts.
10. Performs module-size robustness analysis across top-10, top-15, top-20, top-30, and top-50 definitions.
11. Derives an HBV injury axis from GSE83148 using an ordinal ALT/AST/HBV-DNA injury index.
12. Generates HBV_INJURY top-N and extended gene sets.
13. Exports GSE121248 mixture data for CIBERSORTx.
14. Imports CIBERSORTx cell fractions and fits proliferation + immune-composition adjusted HBV injury models.
15. Computes TCGA-LIHC module scores and Cox proportional-hazards models.
16. Generates manuscript-ready figures and tables.

## Repository structure

```text
R/_shared.R                                      shared functions
config/analysis.yml                              global analysis settings
config/datasets.yml                              cohort registry and label rules
scripts/00_install_packages.R                    install dependencies
scripts/01_dataset_inventory.R                   dataset inventory / inclusion table
scripts/02_download_and_curate_geo.R             GEO download and preprocessing
scripts/02b_audit_gene_symbol_overlap.R          gene-symbol overlap audit
scripts/02c_audit_tissue_labels.R                tissue-label audit
scripts/03_discovery_hallmark_gsva.R             Hallmark GSVA and pathway contrasts
scripts/03b_derive_hepatitis_axes.R              HBV injury axis and top-N gene sets
scripts/04_meta_modules.R                        meta-analysis and module construction
scripts/05_validate_geo_cohorts.R                external GEO validation
scripts/06_tcga_lihc_pipeline.R                  TCGA-LIHC module and Cox models
scripts/07_estimate_adjustment.R                 ESTIMATE exploratory adjustment
scripts/08_cibersortx_export_GSE121248.R         CIBERSORTx mixture export
scripts/09_cibersortx_adjusted_regression_GSE121248.R
                                                   HBV injury adjusted regression
scripts/10_module_size_sensitivity.R             module-size robustness
scripts/11_make_revision_figures.R               Figures 1, 2, 4 and supplement
scripts/11b_make_missing_figures.R               Figures 3 and 5
scripts/12_make_manuscript_tables.R              manuscript tables
scripts/13_summarize_revision_results.R          consolidated result summaries
scripts/run_all.R                                main pipeline launcher
```

## Requirements

R >= 4.3 is recommended.

Install dependencies with:

```bash
Rscript scripts/00_install_packages.R
```

The pipeline uses CRAN and Bioconductor packages including `GEOquery`, `Biobase`, `limma`, `GSVA`, `msigdbr`, `metafor`, `pROC`, `TCGAbiolinks`, `survival`, `broom`, and `tidyverse` packages.

## Quick start

Run the core pipeline:

```bash
Rscript scripts/00_install_packages.R
Rscript scripts/run_all.R
```

`run_all.R` stops before the CIBERSORTx-adjusted regression because CIBERSORTx must be run externally.

For a faster first pass, set the following in `config/analysis.yml`:

```yaml
run_tcga: false
run_estimate: false
```

## CIBERSORTx workflow

First export the mixture file:

```bash
Rscript scripts/08_cibersortx_export_GSE121248.R
```

This creates:

```text
data/external/CIBERSORTx/GSE121248_CIBERSORTx_mixture.txt
```

Run CIBERSORTx externally with:

```text
Mode: Impute Cell Fractions
Signature matrix: LM22
Analysis mode: Relative
Permutations: 100 or 500
Quantile normalization: enabled for microarray input
Batch correction: disabled
Absolute mode: disabled
```

Save the result file as:

```text
data/external/CIBERSORTx/GSE121248_CIBERSORTx_Results.txt
```

Then run:

```bash
Rscript scripts/09_cibersortx_adjusted_regression_GSE121248.R
Rscript scripts/11b_make_missing_figures.R
Rscript scripts/12_make_manuscript_tables.R
Rscript scripts/13_summarize_revision_results.R
```

## Main outputs

Important tables:

```text
results/tables/dataset_inventory.tsv
results/tables/inclusion_exclusion_table.tsv
results/tables/tissue_label_audit.tsv
results/tables/gene_symbol_overlap_audit.tsv
results/tables/discovery_hallmark_combined.tsv
results/tables/meta_analysis_with_heterogeneity.tsv
results/tables/module_gene_panels_top20.tsv
results/tables/geo_validation_summary.tsv
results/tables/module_size_robustness.tsv
results/tables/GSE83148_HBV_INJURY_derivation_full.tsv
results/tables/HBV_INJURY_TOP_2000_gene_set.tsv
results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv
results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv
results/tables/tcga_lihc_survival_cox_models.tsv
```

Manuscript-ready tables are written to:

```text
manuscript/tables/
```

Important figures:

```text
results/figures/Fig1_revised_workflow.pdf/png
results/figures/Fig2_hallmark_dotplot.pdf/png
results/figures/Fig3_discovery_module_boxplots.pdf/png
results/figures/Fig4_multicohort_validation.pdf/png
results/figures/Fig5_HBV_INJURY_adjusted_coefficients.pdf/png
results/figures/SuppFig_module_size_robustness.pdf
```

Copies are also written to:

```text
paper_package/figures/
```

## Manuscript interpretation rules

- `HCCStateScore` is a **bulk transcriptomic tumor-state score**, not an early detection biomarker.
- Tumor/non-tumor AUC values quantify separability in retrospective tissue transcriptomes.
- `HBV_INJURY_TOP_2000` is the compact primary injury representation.
- `HBV_INJURY_EXTENDED_7792` is the extended sensitivity injury representation.
- Persistence of the HBV injury coefficient after E2F/G2M and CIBERSORTx adjustment supports a **residual tumor-associated HBV injury component**, not definitive tumor-cell-intrinsic biology.

## Recommended citation in the manuscript

The repository and generated intermediate files should be cited in the Data Availability statement. Raw data remain available from GEO and TCGA/GDC.
