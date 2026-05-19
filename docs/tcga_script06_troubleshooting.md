# TCGA-LIHC script 06 troubleshooting

## Duplicate `sample_id` in `rownames_to_column()`

`TCGAbiolinks::GDCprepare()` may return `colData(se)` with an existing `sample_id` column. Calling `rownames_to_column("sample_id")` then fails because the requested column name is duplicated. The patched script renames any pre-existing `sample_id` metadata field to `sample_id_coldata` and uses the assay column names as the canonical TCGA sample barcode.

## Missing `sample_type`

Different TCGAbiolinks/GDC releases expose sample annotations under different names (`sample_type`, `definition`, `shortLetterCode`, etc.). The patched script infers tissue from the TCGA barcode sample code first: `01` = tumor and `11` = solid tissue normal, then falls back to text metadata.

## Expected outputs

- `results/tables/tcga_lihc_tissue_label_audit.tsv`
- `results/tables/tcga_lihc_module_scores.tsv`
- `results/tables/tcga_lihc_module_scores_with_clinical.tsv`
- `results/tables/tcga_lihc_clinicopathologic_associations.tsv`
- `results/tables/tcga_lihc_survival_cox_models.tsv`
