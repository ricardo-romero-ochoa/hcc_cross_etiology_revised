source("R/_shared.R")
ensure_dirs(); cfg <- read_config()
steps <- c(
  "scripts/01_dataset_inventory.R",
  "scripts/02_download_and_curate_geo.R",
  "scripts/03_discovery_hallmark_gsva.R",
  "scripts/03b_derive_hepatitis_axes.R",
  "scripts/04_meta_modules.R",
  "scripts/05_validate_geo_cohorts.R",
  "scripts/10_module_size_sensitivity.R",
  "scripts/08_cibersortx_export_GSE121248.R",
  "scripts/11_make_revision_figures.R",
  "scripts/12_make_manuscript_tables.R"
)
if (isTRUE(cfg$run_estimate)) steps <- append(steps, "scripts/07_estimate_adjustment.R", after = 8)
if (isTRUE(cfg$run_tcga)) steps <- append(steps, "scripts/06_tcga_lihc_pipeline.R", after = 8)
for (s in steps) {
  log_msg("Running", s)
  tryCatch(sys.source(s, envir = new.env(parent = globalenv())), error = function(e) {
    writeLines(conditionMessage(e), paste0("results/logs/ERROR_", basename(s), ".log"))
    stop(e)
  })
}
message("Pipeline completed. Run scripts/09_cibersortx_adjusted_regression_GSE121248.R after generating the external CIBERSORTx output file.")
