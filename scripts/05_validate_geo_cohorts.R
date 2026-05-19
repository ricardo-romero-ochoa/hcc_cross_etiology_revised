source("R/_shared.R")
ensure_dirs(); ds <- read_datasets(); panels <- read_module_panels()
validate_one <- function(e) {
  acc <- e$accession
  obj <- read_rds(paste0("data/processed/", acc, "_curated.rds"))
  expr <- obj$expr; meta <- obj$metadata |> mutate(sample_id = as.character(sample_id))
  keep_genes <- c(panels$ProlifHub, panels$HepLoss)
  if (sum(panels$ProlifHub %in% rownames(expr)) < 2 || sum(panels$HepLoss %in% rownames(expr)) < 2) {
    warning("Too few module genes in ", acc); return(NULL)
  }
  scores <- score_modules(expr, panels$ProlifHub, panels$HepLoss)
  readr::write_tsv(scores |> left_join(meta, by = "sample_id"), paste0("results/tables/", acc, "_module_scores.tsv"))
  bind_rows(
    effect_summary(scores, meta, "ProlifHubScore", acc),
    effect_summary(scores, meta, "HepLossScore", acc),
    effect_summary(scores, meta, "HCCStateScore", acc)
  )
}
summary <- map_dfr(ds$validation, validate_one) |> mutate(FDR = p.adjust(p_value, method = "BH"))
readr::write_tsv(summary, "results/tables/geo_validation_summary.tsv")
