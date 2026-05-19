source("R/_shared.R")
suppressPackageStartupMessages({ library(estimate) })
ensure_dirs()
# ESTIMATE expects a GCT file. This script writes gene-symbol expression matrices and runs ESTIMATE where available.
for (acc in c("GSE121248", "GSE41804")) {
  obj <- read_rds(paste0("data/processed/", acc, "_curated.rds"))
  expr <- obj$expr
  gct <- paste0("data/processed/", acc, "_estimate_input.gct")
  filt <- paste0("data/processed/", acc, "_estimate_input_filtered.gct")
  out <- paste0("results/tables/", acc, "_estimate_scores.gct")
  # Write GCT-like matrix. ESTIMATE requires NAME and Description.
  mat <- cbind(NAME = rownames(expr), Description = rownames(expr), as.data.frame(expr, check.names = FALSE))
  con <- file(gct, "w")
  writeLines("#1.2", con)
  writeLines(paste(nrow(expr), ncol(expr), sep = "\t"), con)
  close(con)
  suppressWarnings(readr::write_tsv(as_tibble(mat), gct, append = TRUE))
  tryCatch({
    estimate::filterCommonGenes(input.f = gct, output.f = filt, id = "GeneSymbol")
    estimate::estimateScore(input.ds = filt, output.ds = out, platform = "affymetrix")
    sc <- readr::read_tsv(out, skip = 2, show_col_types = FALSE)
    readr::write_tsv(sc, paste0("results/tables/", acc, "_estimate_scores.tsv"))
  }, error = function(e) {
    writeLines(conditionMessage(e), paste0("results/logs/", acc, "_estimate_error.log"))
  })
}
