source("R/_shared.R")
suppressPackageStartupMessages({ library(msigdbr); library(GSVA); library(limma); library(dplyr) })
ensure_dirs(); cfg <- read_config()

# Robust Hallmark retrieval. Do not use `|> split(x = .$gene_symbol, ...)`:
# the `.` placeholder belongs to magrittr/dplyr pipes, not base R's `|>`.
# This also handles recent msigdbr versions where `category` is deprecated
# in favor of `collection`.
get_hallmark_sets <- function() {
  msig <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  gene_col <- dplyr::case_when(
    "gene_symbol" %in% names(msig) ~ "gene_symbol",
    "db_gene_symbol" %in% names(msig) ~ "db_gene_symbol",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col) || !("gs_name" %in% names(msig))) {
    stop(
      "Could not identify gene-symbol or gene-set columns in msigdbr output. Columns found: ",
      paste(names(msig), collapse = ", ")
    )
  }
  msig <- msig |>
    dplyr::select(gs_name, gene_symbol = dplyr::all_of(gene_col)) |>
    dplyr::mutate(gene_symbol = toupper(gene_symbol)) |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
    dplyr::distinct(gs_name, gene_symbol)
  split(msig$gene_symbol, msig$gs_name)
}

sets <- get_hallmark_sets()
for (acc in c("GSE121248", "GSE41804")) {
  obj <- read_rds(paste0("data/processed/", acc, "_curated.rds"))
  expr <- obj$expr; meta <- obj$metadata
  gs <- run_gsva_scores(expr, sets, method = "gsva", kcdf = "Gaussian",
                        min_size = cfg$gsva$min_size %||% 10,
                        max_size = cfg$gsva$max_size %||% 500,
                        verbose = FALSE)
  fit <- limma_tumor_contrast(gs, meta)
  out <- fit$table |> mutate(dataset = acc, gene_set = gene) |> select(dataset, gene_set, logFC, t, p_value, FDR)
  readr::write_tsv(out, paste0("results/tables/", acc, "_hallmark_limma.tsv"))
  save_rds(list(scores = gs, results = out), paste0("data/processed/", acc, "_hallmark_scores.rds"))
}
bind_rows(lapply(c("GSE121248", "GSE41804"), function(acc) readr::read_tsv(paste0("results/tables/", acc, "_hallmark_limma.tsv"), show_col_types = FALSE))) |>
  readr::write_tsv("results/tables/discovery_hallmark_combined.tsv")
