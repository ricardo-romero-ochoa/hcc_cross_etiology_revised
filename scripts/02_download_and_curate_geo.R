source("R/_shared.R")
ensure_dirs(); cfg <- read_config(); ds <- read_datasets()
entries <- c(ds$discovery, ds$validation)
for (e in entries) {
  acc <- e$accession
  eset <- fetch_geo_eset(acc, e$platform_preferred %||% NULL)
  expr <- collapse_to_symbols(eset)
  pdat <- pData(eset) |> as.data.frame() |> rownames_to_column("sample_id") |> as_tibble()
  tissue <- curate_tissue(pdat |> column_to_rownames("sample_id"), e$tumor_regex, e$nontumor_regex, accession = acc)
  meta <- pdat |> mutate(tissue = tissue, accession = acc, role = e$role, etiology = e$etiology)
  check_tissue_labels(meta, acc)
  # Keep expression columns in the same order as metadata.
  expr <- expr[, meta$sample_id, drop = FALSE]
  save_rds(list(expr = expr, metadata = meta, annotation = annotation(eset)), paste0("data/processed/", acc, "_curated.rds"))
  readr::write_tsv(meta, paste0("results/tables/curated_metadata_", acc, ".tsv"))
  readr::write_tsv(tibble(gene = rownames(expr)), paste0("results/tables/gene_universe_", acc, ".tsv"))
}
# hepatitis-stage cohorts are downloaded and cached for manual axis derivation scripts
for (e in ds$hepatitis_axes) {
  acc <- e$accession
  eset <- fetch_geo_eset(acc, e$platform_preferred %||% NULL)
  expr <- collapse_to_symbols(eset)
  meta <- pData(eset) |> as.data.frame() |> rownames_to_column("sample_id") |> as_tibble() |> mutate(accession = acc, role = e$role, etiology = e$etiology)
  expr <- expr[, meta$sample_id, drop = FALSE]
  save_rds(list(expr = expr, metadata = meta, annotation = annotation(eset)), paste0("data/processed/", acc, "_curated.rds"))
  readr::write_tsv(meta, paste0("results/tables/curated_metadata_", acc, ".tsv"))
}
