# Shared utilities for HCC cross-etiology revision pipeline

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(yaml)
  library(ggplot2)
  library(pROC)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_config <- function(path = "config/analysis.yml") yaml::read_yaml(path)
read_datasets <- function(path = "config/datasets.yml") yaml::read_yaml(path)

ensure_dirs <- function() {
  dirs <- c("data/raw/GEO", "data/processed", "data/external/CIBERSORTx",
            "results/tables", "results/figures", "results/logs")
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  options(timeout = max(600, getOption("timeout", 60)))
}

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  message(msg)
}

geo_series_prefix <- function(accession) {
  sub("[0-9]{1,3}$", "nnn", accession)
}

geo_matrix_dir_url <- function(accession) {
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", geo_series_prefix(accession),
         "/", accession, "/matrix/")
}

list_geo_matrix_urls <- function(accession) {
  base <- geo_matrix_dir_url(accession)
  html <- tryCatch(readLines(base, warn = FALSE), error = function(e) character())
  hits <- unique(unlist(regmatches(html, gregexpr("[A-Za-z0-9_.-]*series_matrix\\.txt\\.gz", html))))
  if (length(hits) == 0) {
    hits <- paste0(accession, "_series_matrix.txt.gz")
  }
  paste0(base, hits)
}

safe_download <- function(url, destfile, retries = 3) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destfile) && file.info(destfile)$size > 0) return(destfile)
  methods <- c("auto", "libcurl", "curl", "wget", "wininet")
  last_error <- NULL
  for (i in seq_len(retries)) {
    for (m in methods) {
      ok <- tryCatch({
        suppressWarnings(utils::download.file(url, destfile = destfile, mode = "wb", method = m, quiet = FALSE))
        file.exists(destfile) && file.info(destfile)$size > 0
      }, error = function(e) { last_error <<- conditionMessage(e); FALSE })
      if (isTRUE(ok)) return(destfile)
      if (file.exists(destfile)) unlink(destfile)
    }
    Sys.sleep(2 * i)
  }
  stop("Could not download ", url, " to ", destfile, ". Last error: ", last_error)
}

fetch_geo_eset <- function(accession, preferred_platform = NULL, destdir = "data/raw/GEO") {
  ensure_dirs()
  destdir <- normalizePath(destdir, mustWork = FALSE)
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Downloading/loading", accession)

  gsets <- tryCatch(
    GEOquery::getGEO(accession, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = destdir),
    error = function(e) {
      log_msg("GEOquery native download failed for", accession, "--", conditionMessage(e))
      NULL
    }
  )

  if (is.null(gsets)) {
    urls <- list_geo_matrix_urls(accession)
    log_msg("Fallback GEO matrix URLs for", accession, ":", paste(basename(urls), collapse = ", "))
    files <- vapply(urls, function(u) safe_download(u, file.path(destdir, basename(u))), character(1))
    gsets <- lapply(files, function(f) GEOquery::getGEO(filename = f, GSEMatrix = TRUE))
  }

  if (!is.list(gsets)) gsets <- list(gsets)
  if (length(gsets) > 1 && !is.null(preferred_platform)) {
    ann <- vapply(gsets, function(x) as.character(annotation(x))[1], character(1))
    idx <- which(ann == preferred_platform | grepl(preferred_platform, basename(names(gsets)), fixed = TRUE))
    if (length(idx) > 0) return(gsets[[idx[1]]])
    log_msg("Preferred platform", preferred_platform, "not found for", accession,
            "; using first matrix with annotation", ann[1])
  }
  gsets[[1]]
}

fetch_gpl_table <- function(gpl_id, destdir = "data/raw/GEO") {
  if (is.null(gpl_id) || is.na(gpl_id) || !grepl("^GPL", gpl_id)) return(NULL)
  log_msg("Fetching platform annotation", gpl_id)
  gpl <- tryCatch(GEOquery::getGEO(gpl_id, AnnotGPL = TRUE, destdir = destdir),
                  error = function(e) {
                    log_msg("Could not fetch", gpl_id, ":", conditionMessage(e)); NULL
                  })
  if (is.null(gpl)) return(NULL)
  as_tibble(GEOquery::Table(gpl))
}

# Robust probe-to-HGNC symbol mapping. GEO platform tables often contain
# several gene-like columns (Gene ID, Entrez Gene, gene_assignment, Gene Symbol).
# Hallmark/GSVA requires HGNC symbols, so avoid numeric Gene ID columns.
clean_symbol <- function(sym_raw) {
  sym <- as.character(sym_raw) |>
    stringr::str_replace("///.*$", "") |>
    stringr::str_replace("//.*$", "") |>
    stringr::str_replace(";.*$", "") |>
    stringr::str_replace(" \\///.*$", "") |>
    stringr::str_replace(" \\(.*$", "") |>
    stringr::str_trim()
  sym[sym == "" | sym == "---" | sym == "NA" | is.na(sym)] <- NA_character_
  sym
}

symbol_like_fraction <- function(x) {
  sym <- clean_symbol(x)
  sym <- sym[!is.na(sym)]
  if (length(sym) == 0) return(0)
  mean(stringr::str_detect(sym, "[A-Za-z]") & !stringr::str_detect(sym, "^[0-9.]+$"))
}

find_symbol_column <- function(df) {
  nm <- colnames(df)
  preferred_patterns <- c(
    "^Gene[._ ]?Symbol$", "^GENE[._ ]?SYMBOL$", "^gene[._ ]?symbol$",
    "^Symbol$", "^SYMBOL$", "^symbol$", "^ILMN_Gene$"
  )
  candidates <- unique(unlist(lapply(preferred_patterns, function(p) grep(p, nm, ignore.case = FALSE, value = TRUE))))
  if (length(candidates) == 0) candidates <- grep("symbol", nm, ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0) candidates <- grep("gene_assignment|assignment|gene.?symbol", nm, ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0) return(NA_character_)
  bad <- grep("entrez|gene.?id|identifier|gb_acc|refseq|unigene", candidates, ignore.case = TRUE, value = TRUE)
  candidates <- setdiff(candidates, bad)
  if (length(candidates) == 0) return(NA_character_)
  frac <- vapply(candidates, function(cc) symbol_like_fraction(df[[cc]]), numeric(1))
  candidates <- candidates[order(frac, decreasing = TRUE)]
  if (length(candidates) == 0 || max(frac, na.rm = TRUE) < 0.25) return(NA_character_)
  candidates[1]
}

safe_symbol_map <- function(eset) {
  fdat <- fData(eset)
  if (nrow(fdat) > 0) {
    hit <- find_symbol_column(fdat)
    if (!is.na(hit)) {
      sym <- clean_symbol(fdat[[hit]])
      return(tibble(probe_id = rownames(fdat), symbol = sym))
    }
  }

  gpl_id <- as.character(annotation(eset))[1]
  gpl_tab <- fetch_gpl_table(gpl_id)
  if (!is.null(gpl_tab)) {
    id_col <- intersect(c("ID", "ID_REF", "Probe Set ID", "ProbeID", "PROBE_ID"), colnames(gpl_tab))
    if (length(id_col) == 0) id_col <- colnames(gpl_tab)[1]
    sym_col <- find_symbol_column(gpl_tab)
    if (!is.na(sym_col)) {
      map <- tibble(probe_id = as.character(gpl_tab[[id_col[1]]]), symbol = clean_symbol(gpl_tab[[sym_col]]))
      return(tibble(probe_id = rownames(eset), symbol = map$symbol[match(rownames(eset), map$probe_id)]))
    }
  }
  stop("No gene-symbol-like annotation column was found for ", gpl_id,
       ". Check platform annotation or provide a custom probe-to-symbol map.")
}

collapse_to_symbols <- function(eset) {
  expr <- exprs(eset)
  map <- safe_symbol_map(eset)
  keep <- !is.na(map$symbol) & map$symbol != ""
  if (sum(keep) < 1000) {
    stop("Probe-to-symbol mapping produced too few valid HGNC symbols for platform ",
         as.character(annotation(eset))[1], ". Check fData/GPL annotation.")
  }
  expr <- expr[map$probe_id[keep], , drop = FALSE]
  map <- map[keep, ]
  map$symbol <- toupper(map$symbol)
  avg <- rowMeans(expr, na.rm = TRUE)
  o <- order(map$symbol, -avg)
  expr_o <- expr[o, , drop = FALSE]
  map_o <- map[o, ]
  keep_first <- !duplicated(map_o$symbol)
  out <- expr_o[keep_first, , drop = FALSE]
  rownames(out) <- map_o$symbol[keep_first]
  attr(out, "n_mapped_probes") <- nrow(expr)
  attr(out, "n_unique_symbols") <- nrow(out)
  out
}

# Build a compact sample-level text field for phenotype parsing. Avoid columns that
# are constant across all samples (e.g., series title, study-wide HCC description),
# because they can make every sample look like "tumor".
metadata_text <- function(pdat) {
  pdat <- as.data.frame(pdat)
  preferred <- grep(
    "^(title|source_name_ch1|characteristics_ch1|description|sample_type|tissue|disease|group|condition|status)",
    colnames(pdat), ignore.case = TRUE, value = TRUE
  )
  if (length(preferred) == 0) preferred <- colnames(pdat)
  x <- pdat[, preferred, drop = FALSE]

  nunique <- vapply(x, function(z) length(unique(stats::na.omit(as.character(z)))), integer(1))
  keep_cols <- names(nunique)[nunique > 1]
  if (length(keep_cols) > 0) x <- x[, keep_cols, drop = FALSE]

  apply(x, 1, function(r) {
    r <- as.character(r)
    r <- r[!is.na(r) & trimws(r) != ""]
    paste(r, collapse = " | ")
  })
}

# Tumor/non-tumor label parser. Non-tumor is deliberately evaluated first because
# strings such as "non-tumor" contain the token "tumor". Tumor-first parsing can
# label every sample as tumor and then limma reports: "Coefficients not estimable:
# tissuetumor".
curate_tissue <- function(pdat, tumor_regex, nontumor_regex, accession = NA_character_) {
  txt <- metadata_text(pdat)
  txt_l <- stringr::str_to_lower(txt)
  tissue <- rep(NA_character_, length(txt_l))

  strong_non <- paste(c(
    "non[-_ ]?tumou?r", "non[-_ ]?cancer", "noncancerous", "non[-_ ]?neoplastic",
    "adjacent", "para[-_ ]?tumou?r", "peri[-_ ]?tumou?r", "surrounding",
    "paired normal", "normal liver", "normal tissue", "tumou?r[-_ ]?free",
    "cirrhotic", "cirrhosis", "healthy"
  ), collapse = "|")

  strong_tumor <- paste(c(
    "primary tumou?r", "tumou?r tissue", "hcc tumou?r", "hepatocellular carcinoma",
    "cancerous tissue", "carcinoma tissue", "\\btumou?r\\b", "\\bcancer\\b"
  ), collapse = "|")

  is_non <- stringr::str_detect(txt_l, strong_non) | stringr::str_detect(txt, nontumor_regex)
  tissue[is_non] <- "non_tumor"

  is_tumor <- (stringr::str_detect(txt_l, strong_tumor) | stringr::str_detect(txt, tumor_regex)) & !is_non
  tissue[is_tumor] <- "tumor"

  if (!is.na(accession)) {
    acc <- as.character(accession)
    if (acc %in% c("GSE121248", "GSE41804", "GSE14520", "GSE25097", "GSE76427", "GSE36376", "GSE57957", "GSE45267", "GSE112790")) {
      unknown <- is.na(tissue)
      tissue[unknown & stringr::str_detect(txt_l[unknown], "\\bt\\b|tumou?r|carcinoma|cancer")] <- "tumor"
      unknown <- is.na(tissue)
      tissue[unknown & stringr::str_detect(txt_l[unknown], "\\bn\\b|normal|adjacent|non|cirrh|healthy")] <- "non_tumor"
    }
  }

  tissue
}

check_tissue_labels <- function(metadata, accession = NA_character_) {
  tab <- table(metadata$tissue, useNA = "ifany")
  msg <- paste(names(tab), as.integer(tab), sep = "=", collapse = ", ")
  message("[labels] ", accession, ": ", msg)
  if (!all(c("tumor", "non_tumor") %in% metadata$tissue)) {
    preview_cols <- intersect(c("sample_id", "title", "source_name_ch1", grep("^characteristics_ch1", names(metadata), value = TRUE)), names(metadata))
    preview <- utils::capture.output(print(utils::head(metadata[, preview_cols, drop = FALSE], 12)))
    stop(
      "Tissue parsing failed for ", accession, ": both tumor and non_tumor are required. Counts: ", msg,
      "\nInspect results/tables/curated_metadata_", accession, ".tsv or adjust config/datasets.yml regex.\n",
      paste(preview, collapse = "\n")
    )
  }
  invisible(tab)
}

make_design_tissue <- function(metadata) {
  metadata$tissue <- factor(metadata$tissue, levels = c("non_tumor", "tumor"))
  model.matrix(~ tissue, data = metadata)
}

limma_tumor_contrast <- function(expr, metadata) {
  keep <- !is.na(metadata$tissue) & metadata$tissue %in% c("tumor", "non_tumor")
  expr <- expr[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  if (!all(c("tumor", "non_tumor") %in% metadata$tissue)) {
    tab <- table(metadata$tissue, useNA = "ifany")
    stop("Cannot fit tumor contrast: tissue labels do not contain both tumor and non_tumor. Counts: ",
         paste(names(tab), as.integer(tab), sep = "=", collapse = ", "))
  }
  design <- make_design_tissue(metadata)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "tissuetumor", number = Inf, sort.by = "none") |>
    rownames_to_column("gene") |>
    as_tibble() |>
    rename(logFC = logFC, p_value = P.Value, FDR = adj.P.Val, t = t)
  list(table = tt, fit = fit, expr = expr, metadata = metadata)
}

zscore_rows <- function(expr) {
  t(scale(t(expr)))
}

score_module <- function(expr, genes, min_genes = 2) {
  genes <- toupper(unique(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]
  rownames(expr) <- toupper(rownames(expr))
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < min_genes) {
    warning("Only ", length(genes), " module genes overlap the expression matrix; minimum required is ", min_genes, ". Returning NA scores.")
    return(rep(NA_real_, ncol(expr)))
  }
  z <- zscore_rows(expr[genes, , drop = FALSE])
  as.numeric(colMeans(z, na.rm = TRUE))
}

score_modules <- function(expr, prolif_genes, heploss_genes) {
  pro <- score_module(expr, prolif_genes)
  hep <- score_module(expr, heploss_genes)
  tibble(sample_id = colnames(expr),
         ProlifHubScore = as.numeric(pro),
         HepLossScore = as.numeric(hep),
         HCCStateScore = ProlifHubScore - HepLossScore)
}



# Version-safe GSVA wrapper.
# Bioconductor GSVA changed its API: recent versions require a *Param object
# and no longer dispatch gsva() on a plain matrix. Older versions use
# gsva(expr, gene_sets, method=..., kcdf=...). This wrapper supports both APIs.
filter_gene_sets <- function(gene_sets, expr, min_size = 10, max_size = 500) {
  rownames(expr) <- toupper(rownames(expr))
  gene_sets <- lapply(gene_sets, function(g) intersect(toupper(unique(g)), rownames(expr)))
  sizes <- lengths(gene_sets)
  gene_sets[sizes >= min_size & sizes <= max_size]
}

gene_set_overlap_diagnostics <- function(expr, gene_sets) {
  expr_genes <- toupper(rownames(expr))
  gs_genes <- unique(toupper(unlist(gene_sets, use.names = FALSE)))
  n_overlap <- length(intersect(expr_genes, gs_genes))
  probe_sample <- expr_genes[seq_len(min(1000, length(expr_genes)))]
  looks_like_probes <- mean(stringr::str_detect(probe_sample, "(_AT$|_S_AT$|_X_AT$|^[0-9]+$|^ILMN_|^ENSG)"), na.rm = TRUE)
  list(
    n_expr_genes = length(expr_genes),
    n_gene_set_genes = length(gs_genes),
    n_overlap = n_overlap,
    first_expr_ids = paste(utils::head(expr_genes, 10), collapse = ", "),
    looks_like_probes = looks_like_probes
  )
}

run_gsva_scores <- function(expr, gene_sets, method = "gsva", kcdf = "Gaussian",
                            min_size = 10, max_size = 500, verbose = FALSE) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("The GSVA package is required. Install it with BiocManager::install('GSVA').")
  }
  expr <- as.matrix(expr)
  rownames(expr) <- toupper(rownames(expr))
  storage.mode(expr) <- "numeric"
  diag0 <- gene_set_overlap_diagnostics(expr, gene_sets)
  gene_sets <- filter_gene_sets(gene_sets, expr, min_size = min_size, max_size = max_size)
  if (length(gene_sets) == 0) {
    stop(
      "No gene sets remained after filtering to genes present in the expression matrix. ",
      "This almost always means that the processed expression matrix still uses probe/Entrez/Ensembl IDs ",
      "instead of HGNC symbols, or that script 02 was not rerun after the annotation patch. ",
      "Overlap before size filtering: ", diag0$n_overlap, " / ", diag0$n_gene_set_genes, "; ",
      "expression rows: ", diag0$n_expr_genes, "; first expression row IDs: ", diag0$first_expr_ids, "; ",
      "probe-like fraction among first rows: ", signif(diag0$looks_like_probes, 3), ". ",
      "Delete data/processed/*_curated.rds and rerun scripts/02_download_and_curate_geo.R."
    )
  }

  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    ctor <- get("gsvaParam", envir = asNamespace("GSVA"))
    fml <- names(formals(ctor))
    args <- list()
    if ("exprData" %in% fml) args$exprData <- expr else if ("expr" %in% fml) args$expr <- expr else args[[1]] <- expr
    if ("geneSets" %in% fml) args$geneSets <- gene_sets else if ("gset.idx.list" %in% fml) args$gset.idx.list <- gene_sets else args[[length(args) + 1]] <- gene_sets
    if ("kcdf" %in% fml) args$kcdf <- kcdf
    if ("minSize" %in% fml) args$minSize <- min_size
    if ("maxSize" %in% fml) args$maxSize <- max_size
    param <- tryCatch(do.call(ctor, args), error = function(e) e)
    if (!inherits(param, "error")) {
      res <- tryCatch(GSVA::gsva(param, verbose = verbose), error = function(e1) tryCatch(GSVA::gsva(param), error = function(e2) e2))
      if (!inherits(res, "error")) return(as.matrix(res))
      message("New GSVA API call failed: ", conditionMessage(res))
    } else {
      message("gsvaParam construction failed: ", conditionMessage(param))
    }
  }

  old_call <- tryCatch(
    GSVA::gsva(expr, gene_sets, method = method, kcdf = kcdf, min.sz = min_size, max.sz = max_size, verbose = verbose),
    error = function(e) e
  )
  if (!inherits(old_call, "error")) return(as.matrix(old_call))

  stop("GSVA scoring failed under both the new Param API and the old matrix API. Installed GSVA version: ",
       as.character(utils::packageVersion("GSVA")), ". Last error: ", conditionMessage(old_call))
}


effect_summary <- function(scores, metadata, score_col, dataset) {
  dat <- scores |> left_join(metadata |> select(sample_id, tissue), by = "sample_id") |>
    filter(tissue %in% c("tumor", "non_tumor"), !is.na(.data[[score_col]]))
  x <- dat |> filter(tissue == "tumor") |> pull(all_of(score_col))
  y <- dat |> filter(tissue == "non_tumor") |> pull(all_of(score_col))
  delta <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  sp <- sqrt(((length(x)-1)*var(x) + (length(y)-1)*var(y))/(length(x)+length(y)-2))
  d <- delta / sp
  p <- t.test(x, y)$p.value
  auc <- NA_real_; auc_low <- NA_real_; auc_high <- NA_real_
  if (length(unique(dat$tissue)) == 2) {
    rr <- pROC::roc(response = dat$tissue, predictor = dat[[score_col]],
                    levels = c("non_tumor", "tumor"), quiet = TRUE, direction = "<")
    auc <- as.numeric(pROC::auc(rr))
    ci <- suppressMessages(pROC::ci.auc(rr))
    auc_low <- as.numeric(ci[1]); auc_high <- as.numeric(ci[3])
  }
  tibble(dataset = dataset, score = score_col, n_tumor = length(x), n_non_tumor = length(y),
         delta_tumor_minus_nontumor = delta, cohens_d = d, p_value = p,
         AUC = auc, AUC_low = auc_low, AUC_high = auc_high)
}

save_rds <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); saveRDS(x, path); invisible(x) }
read_rds <- function(path) readRDS(path)

read_module_panels <- function(path = "results/tables/module_gene_panels_top20.tsv") {
  tab <- readr::read_tsv(path, show_col_types = FALSE)
  list(ProlifHub = tab |> filter(module == "ProlifHub") |> pull(gene),
       HepLoss = tab |> filter(module == "HepLoss") |> pull(gene))
}

clr_transform <- function(mat, pseudo = 1e-6) {
  mat <- as.matrix(mat)
  mat <- mat + pseudo
  logmat <- log(mat)
  sweep(logmat, 1, rowMeans(logmat), FUN = "-")
}

extract_tissue_coef <- function(fit_obj, term = "tissuetumor") {
  sm <- summary(fit_obj)$coefficients
  if (!term %in% rownames(sm)) return(tibble(term = term, estimate = NA_real_, std_error = NA_real_, t = NA_real_, p_value = NA_real_))
  tibble(term = term, estimate = sm[term, "Estimate"], std_error = sm[term, "Std. Error"],
         t = sm[term, "t value"], p_value = sm[term, "Pr(>|t|)"])
}
