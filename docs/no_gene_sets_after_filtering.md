# Fix: “No gene sets remained after filtering”

This error means that the expression matrix row names do not overlap with MSigDB/Hallmark HGNC gene symbols. The usual cause is an incorrect GEO platform annotation column, for example using `Gene ID` / Entrez IDs or Affymetrix probe IDs instead of `Gene Symbol`.

## Required cleanup after this patch

Delete old processed matrices because they may have been saved with the wrong row identifiers:

```bash
rm -f data/processed/*_curated.rds
rm -f results/tables/gene_universe_*.tsv
```

Then rerun:

```bash
Rscript scripts/02_download_and_curate_geo.R
Rscript scripts/02b_audit_gene_symbol_overlap.R
Rscript scripts/03_discovery_hallmark_gsva.R
```

The audit writes:

```text
results/tables/gene_symbol_overlap_audit.tsv
```

For GPL570 cohorts, the overlap with Hallmark genes should be large, not zero. If the `first_expression_ids` column contains values like `1007_s_at` or purely numeric IDs, the matrix was not collapsed to HGNC symbols and script 02 must be rerun.
