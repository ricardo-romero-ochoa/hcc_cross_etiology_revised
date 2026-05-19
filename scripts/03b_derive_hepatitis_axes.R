# ============================================================
# 03b_derive_hepatitis_axes.R
# Derive HBV_INJURY programs from GSE83148 using an ordinal
# ALT/AST/HBV-DNA injury index.
#
# Outputs:
#   results/tables/GSE83148_metadata_audit.tsv
#   results/tables/GSE83148_HBV_INJURY_raw_components.tsv
#   results/tables/GSE83148_HBV_INJURY_index_components.tsv
#   results/tables/GSE83148_HBV_INJURY_derivation_full.tsv
#   results/tables/HBV_INJURY_TOP_200_gene_set.tsv
#   results/tables/HBV_INJURY_TOP_500_gene_set.tsv
#   results/tables/HBV_INJURY_TOP_1000_gene_set.tsv
#   results/tables/HBV_INJURY_TOP_2000_gene_set.tsv
#   results/tables/HBV_INJURY_TOP_5000_gene_set.tsv
#   results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv
#   results/tables/HBV_INJURY_topN_extended_summary.tsv
#   results/tables/HBV_INJURY_gene_set.tsv   backward-compatible extended set
# ============================================================

source("R/_shared.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(limma)
  library(stringr)
  library(purrr)
})

ensure_dirs()

obj_path <- "data/processed/GSE83148_curated.rds"
if (!file.exists(obj_path)) {
  stop("Missing ", obj_path, "\nRun scripts/02_download_and_curate_geo.R first and ensure GSE83148 is in config/datasets.yml.")
}

obj <- read_rds(obj_path)
expr <- obj$expr
meta <- obj$metadata |> mutate(sample_id = as.character(sample_id))

# Metadata audit for transparency and reviewer response.
meta_audit <- tibble(
  column = colnames(meta),
  n_non_missing = vapply(meta, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1)),
  n_unique = vapply(meta, function(x) length(unique(na.omit(as.character(x)))), integer(1)),
  example_values = vapply(meta, function(x) {
    vals <- unique(na.omit(as.character(x)))
    vals <- vals[vals != ""]
    paste(head(vals, 8), collapse = " | ")
  }, character(1))
)
readr::write_tsv(meta_audit, "results/tables/GSE83148_metadata_audit.tsv")

find_one_col <- function(meta, patterns, label) {
  hits <- unique(unlist(lapply(patterns, function(p) grep(p, colnames(meta), value = TRUE, ignore.case = TRUE))))
  if (length(hits) == 0) {
    stop("Could not find column for ", label, ". Inspect results/tables/GSE83148_metadata_audit.tsv.")
  }
  hits[1]
}

alt_col <- find_one_col(meta, c("^alt", "alt:ch1", "alanine"), "ALT")
ast_col <- find_one_col(meta, c("^ast", "ast:ch1", "aspartate"), "AST")
dna_col <- find_one_col(meta, c("hbv.*dna", "hbv-dna", "viral.*load"), "HBV-DNA")

message("[03b] ALT column: ", alt_col)
message("[03b] AST column: ", ast_col)
message("[03b] HBV-DNA column: ", dna_col)

# GSE83148-specific ordinal coding:
# NON = 0, <= lower cutoff = 1, > higher cutoff = 2.
recode_threshold_ordinal <- function(x) {
  xs <- str_to_upper(str_squish(as.character(x)))
  dplyr::case_when(
    is.na(xs) | xs == "" ~ NA_real_,
    xs == "NON" ~ 0,
    str_detect(xs, "^<=") ~ 1,
    str_detect(xs, "^<") ~ 1,
    str_detect(xs, "^>=") ~ 2,
    str_detect(xs, "^>") ~ 2,
    TRUE ~ NA_real_
  )
}

component_df <- meta |>
  transmute(
    sample_id,
    ALT_ord = recode_threshold_ordinal(.data[[alt_col]]),
    AST_ord = recode_threshold_ordinal(.data[[ast_col]]),
    HBVDNA_ord = recode_threshold_ordinal(.data[[dna_col]])
  )

readr::write_tsv(component_df, "results/tables/GSE83148_HBV_INJURY_raw_components.tsv")

z_safe <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

injury_index <- component_df |>
  mutate(
    ALT_z = z_safe(ALT_ord),
    AST_z = z_safe(AST_ord),
    HBVDNA_z = z_safe(HBVDNA_ord),
    HBV_INJURY_INDEX = rowMeans(across(c(ALT_z, AST_z, HBVDNA_z)), na.rm = TRUE),
    HBV_INJURY_INDEX = ifelse(is.nan(HBV_INJURY_INDEX), NA_real_, HBV_INJURY_INDEX)
  ) |>
  select(sample_id, ALT_ord, AST_ord, HBVDNA_ord, ALT_z, AST_z, HBVDNA_z, HBV_INJURY_INDEX)

readr::write_tsv(injury_index, "results/tables/GSE83148_HBV_INJURY_index_components.tsv")

common <- intersect(colnames(expr), injury_index$sample_id)
expr2 <- expr[, common, drop = FALSE]
inj <- injury_index |>
  filter(sample_id %in% common) |>
  arrange(match(sample_id, common))
stopifnot(all(colnames(expr2) == inj$sample_id))

usable <- is.finite(inj$HBV_INJURY_INDEX)
expr2 <- expr2[, usable, drop = FALSE]
inj <- inj[usable, , drop = FALSE]

if (ncol(expr2) < 6) stop("Too few samples with usable HBV_INJURY_INDEX.")
message("[03b] Samples with usable HBV_INJURY_INDEX: ", ncol(expr2))

design <- model.matrix(~ HBV_INJURY_INDEX, data = inj)
fit <- limma::lmFit(expr2, design)
fit <- limma::eBayes(fit)

tt_raw <- limma::topTable(fit, coef = "HBV_INJURY_INDEX", number = Inf, sort.by = "P") |>
  rownames_to_column("gene") |>
  as_tibble()

effect_col <- dplyr::case_when(
  "logFC" %in% colnames(tt_raw) ~ "logFC",
  "HBV_INJURY_INDEX" %in% colnames(tt_raw) ~ "HBV_INJURY_INDEX",
  TRUE ~ NA_character_
)
if (is.na(effect_col)) stop("Could not identify effect column in topTable output. Columns: ", paste(colnames(tt_raw), collapse = ", "))

tt <- tt_raw |>
  mutate(
    beta_injury_index = .data[[effect_col]],
    p_value = .data[["P.Value"]],
    FDR = .data[["adj.P.Val"]]
  ) |>
  select(gene, beta_injury_index, AveExpr, t, p_value, FDR, B, everything())

readr::write_tsv(tt, "results/tables/GSE83148_HBV_INJURY_derivation_full.tsv")

ranked <- tt |>
  filter(beta_injury_index > 0, FDR < 0.10) |>
  arrange(FDR, desc(abs(t)), desc(beta_injury_index)) |>
  mutate(rank = row_number(), selection_rule = "positive_beta_FDR_lt_0.10_ranked_by_FDR_and_moderated_t")

if (nrow(ranked) < 10) stop("Too few positive FDR<0.10 injury-associated genes: ", nrow(ranked))

# Extended set; the historical revised analysis had 7,792 genes. The compatibility
# filename keeps the manuscript pipeline stable even if reruns produce a slightly
# different count due to annotation/package changes.
extended <- ranked |>
  mutate(injury_set = paste0("HBV_INJURY_EXTENDED_", nrow(ranked))) |>
  select(injury_set, rank, gene, beta_injury_index, t, p_value, FDR, selection_rule)
readr::write_tsv(extended, paste0("results/tables/HBV_INJURY_EXTENDED_", nrow(ranked), "_gene_set.tsv"))
readr::write_tsv(extended, "results/tables/HBV_INJURY_EXTENDED_7792_gene_set.tsv")
readr::write_tsv(extended, "results/tables/HBV_INJURY_gene_set.tsv")

for (n in c(200, 500, 1000, 2000, 5000)) {
  out <- ranked |>
    slice_head(n = min(n, nrow(ranked))) |>
    mutate(
      injury_set = paste0("HBV_INJURY_TOP_", n),
      selection_rule = paste0("top_", n, "_positive_beta_FDR_lt_0.10_ranked")
    ) |>
    select(injury_set, rank, gene, beta_injury_index, t, p_value, FDR, selection_rule)
  readr::write_tsv(out, paste0("results/tables/HBV_INJURY_TOP_", n, "_gene_set.tsv"))
}

summary_tbl <- bind_rows(
  purrr::map_dfr(c(200, 500, 1000, 2000, 5000), function(n) {
    readr::read_tsv(paste0("results/tables/HBV_INJURY_TOP_", n, "_gene_set.tsv"), show_col_types = FALSE)
  }),
  extended
) |>
  group_by(injury_set, selection_rule) |>
  summarise(
    n_genes = n(),
    min_FDR = min(FDR, na.rm = TRUE),
    max_FDR = max(FDR, na.rm = TRUE),
    median_beta = median(beta_injury_index, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_tsv(summary_tbl, "results/tables/HBV_INJURY_topN_extended_summary.tsv")

message("[03b] Wrote HBV_INJURY derivation and top-N/extended gene sets. Extended n = ", nrow(extended))
