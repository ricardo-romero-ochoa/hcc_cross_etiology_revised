# `meta_p_FE` not found in script 04

This error occurs when the meta-analysis helper returns an empty tibble or rows without the expected fixed-effect p-value column. In practice this can happen when all genes fail the standard-error checks, usually because one of the discovery contrasts produced missing or zero moderated t-statistics, or because no valid common genes remained after preprocessing.

The patched `scripts/04_meta_modules.R` always returns a typed one-row result for every gene, even when a gene is not estimable. It adds a `meta_status` column with values:

- `ok`: fixed-effect meta-analysis was estimated;
- `not_estimable`: missing, zero, non-finite, or invalid statistics prevented meta-analysis.

Useful files to inspect after running script 04:

- `results/tables/discovery_gene_limma_wide.tsv`
- `results/tables/meta_analysis_with_heterogeneity.tsv`

If `meta_status == ok` is rare or absent, check tissue labels and gene-symbol mapping in the curated RDS objects.
