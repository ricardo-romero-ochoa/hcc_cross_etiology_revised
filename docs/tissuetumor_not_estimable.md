# `Coefficients not estimable: tissuetumor`

This warning means the model matrix contains only one usable tissue class, or `tissuetumor` is collinear with the intercept. In this pipeline the usual cause is phenotype parsing: non-tumor samples were matched by the tumor regex because strings such as `non-tumor` contain `tumor`, or because a constant study-level HCC description was used for every sample.

The v7 patch fixes this by:

1. Building phenotype text only from sample-varying columns.
2. Prioritising non-tumor labels before tumor labels.
3. Adding `check_tissue_labels()` after curation.
4. Stopping before limma when both classes are not present.

After updating, delete stale curated files and rerun:

```bash
rm -f data/processed/*_curated.rds
rm -f results/tables/curated_metadata_*.tsv
Rscript scripts/02_download_and_curate_geo.R
Rscript scripts/02c_audit_tissue_labels.R
Rscript scripts/03_discovery_hallmark_gsva.R
Rscript scripts/04_meta_modules.R
```

Inspect:

```text
results/tables/tissue_label_audit.tsv
results/tables/curated_metadata_GSE121248.tsv
results/tables/curated_metadata_GSE41804.tsv
```

Both discovery cohorts must contain non-zero `tumor` and `non_tumor` counts before scripts 03 and 04 can work.
