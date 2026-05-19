source("R/_shared.R")
ensure_dirs(); ds <- read_datasets()
entries <- c(ds$discovery, ds$validation)
out <- lapply(entries, function(e) {
  acc <- e$accession
  f <- paste0("data/processed/", acc, "_curated.rds")
  if (!file.exists(f)) {
    return(tibble::tibble(accession = acc, status = "missing_curated_rds", tissue = NA_character_, n = NA_integer_))
  }
  obj <- readr::read_rds(f)
  tab <- as.data.frame(table(obj$metadata$tissue, useNA = "ifany"))
  names(tab) <- c("tissue", "n")
  tibble::as_tibble(tab) |> dplyr::mutate(accession = acc, status = "ok", .before = 1)
}) |> dplyr::bind_rows()
readr::write_tsv(out, "results/tables/tissue_label_audit.tsv")
print(out)
