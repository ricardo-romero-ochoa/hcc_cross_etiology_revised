source("R/_shared.R")
suppressPackageStartupMessages({ library(msigdbr); library(dplyr); library(readr); library(tibble) })
ensure_dirs()

get_hallmark_sets_for_audit <- function() {
  msig <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  gene_col <- if ("gene_symbol" %in% names(msig)) "gene_symbol" else if ("db_gene_symbol" %in% names(msig)) "db_gene_symbol" else stop("No gene symbol column in msigdbr output")
  msig <- msig |>
    select(gs_name, gene_symbol = all_of(gene_col)) |>
    mutate(gene_symbol = toupper(gene_symbol)) |>
    filter(!is.na(gene_symbol), gene_symbol != "") |>
    distinct(gs_name, gene_symbol)
  split(msig$gene_symbol, msig$gs_name)
}

sets <- get_hallmark_sets_for_audit()
files <- list.files("data/processed", pattern = "_curated\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No curated RDS files found. Run scripts/02_download_and_curate_geo.R first.")
res <- lapply(files, function(f) {
  obj <- readRDS(f)
  expr <- obj$expr
  d <- gene_set_overlap_diagnostics(expr, sets)
  tibble(
    file = basename(f),
    accession = sub("_curated\\.rds$", "", basename(f)),
    n_expression_rows = d$n_expr_genes,
    n_hallmark_unique_genes = d$n_gene_set_genes,
    n_overlap = d$n_overlap,
    first_expression_ids = d$first_expr_ids,
    probe_like_fraction_first_1000 = d$looks_like_probes
  )
}) |> bind_rows()
write_tsv(res, "results/tables/gene_symbol_overlap_audit.tsv")
print(res)
