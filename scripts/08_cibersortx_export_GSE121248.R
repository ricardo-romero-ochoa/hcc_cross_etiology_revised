source("R/_shared.R")
ensure_dirs(); cfg <- read_config()
obj <- read_rds("data/processed/GSE121248_curated.rds")
expr <- obj$expr
# CIBERSORTx mixture format: first column is gene symbol, subsequent columns are samples.
# For Affymetrix processed matrices, the values are usually log2/RMA. CIBERSORTx accepts normalized microarray input;
# use quantile normalization in CIBERSORTx unless using RNA-seq.
mixture <- as.data.frame(expr, check.names = FALSE) |> rownames_to_column("GeneSymbol")
readr::write_tsv(mixture, cfg$cibersortx$mixture_out)
readr::write_tsv(obj$metadata, "data/external/CIBERSORTx/GSE121248_CIBERSORTx_sample_metadata.tsv")
message("Wrote CIBERSORTx mixture file: ", cfg$cibersortx$mixture_out)
message("After running CIBERSORTx, save output as: ", cfg$cibersortx$results_in)
