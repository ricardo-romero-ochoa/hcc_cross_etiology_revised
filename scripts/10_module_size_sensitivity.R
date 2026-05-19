source("R/_shared.R")
ensure_dirs(); cfg <- read_config(); ds <- read_datasets()
validate_size <- function(n, e) {
  panel_path <- paste0("results/tables/module_gene_panels_top", n, ".tsv")
  if (!file.exists(panel_path)) return(NULL)
  panels <- read_module_panels(panel_path)
  acc <- e$accession
  obj <- read_rds(paste0("data/processed/", acc, "_curated.rds"))
  if (sum(panels$ProlifHub %in% rownames(obj$expr)) < 2 || sum(panels$HepLoss %in% rownames(obj$expr)) < 2) return(NULL)
  scores <- score_modules(obj$expr, panels$ProlifHub, panels$HepLoss)
  effect_summary(scores, obj$metadata, "HCCStateScore", acc) |> mutate(module_size = n)
}
rob <- map_dfr(cfg$module_sizes, function(n) map_dfr(c(ds$discovery, ds$validation), ~validate_size(n, .x)))
readr::write_tsv(rob, "results/tables/module_size_robustness.tsv")
