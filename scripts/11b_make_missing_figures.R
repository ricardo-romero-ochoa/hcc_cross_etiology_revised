# ============================================================
# 11b_make_missing_figures.R
# Missing revised figures for HCC manuscript
#
# Generates:
#   paper_package/figures/Fig5_HBV_INJURY_adjusted_coefficients.pdf/png
#   paper_package/figures/Fig3_discovery_module_boxplots.pdf/png
#
# Required upstream outputs:
#   results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv
#   data/processed/GSE121248_curated.rds
#   data/processed/GSE41804_curated.rds
#   results/tables/module_gene_panels_top20.tsv
#
# Optional fallback files:
#   results/tables/gse121248_cibersortx_adjusted_regression.tsv
#   results/modules/prolifhub_genes_top20.txt
#   results/modules/heploss_genes_top20.txt
# ============================================================

source("R/_shared.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(broom)
  library(rlang)
})

# Use whatever directory helper exists in the local repo.
if (exists("ensure_dirs")) {
  ensure_dirs()
} else if (exists("make_dirs")) {
  make_dirs()
} else {
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  dir.create("paper_package/figures", recursive = TRUE, showWarnings = FALSE)
}

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_package/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

p_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 1e-4 ~ "p < 1e-4",
    p < 0.001 ~ paste0("p = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("p = ", signif(p, 3))
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

safe_cohens_d <- function(x_tumor, x_non) {
  x_tumor <- x_tumor[is.finite(x_tumor)]
  x_non <- x_non[is.finite(x_non)]
  if (length(x_tumor) < 2 || length(x_non) < 2) return(NA_real_)
  s1 <- stats::sd(x_tumor)
  s0 <- stats::sd(x_non)
  sp <- sqrt(((length(x_tumor) - 1) * s1^2 + (length(x_non) - 1) * s0^2) /
               (length(x_tumor) + length(x_non) - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x_tumor) - mean(x_non)) / sp
}

score_module_local <- function(expr_gene, genes, min_genes = 5) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  g <- intersect(genes, rownames(expr_gene))
  if (length(g) < min_genes) {
    warning("Only ", length(g), " module genes overlap expression matrix. Returning NA.")
    return(rep(NA_real_, ncol(expr_gene)))
  }
  z <- t(scale(t(expr_gene[g, , drop = FALSE])))
  as.numeric(colMeans(z, na.rm = TRUE))
}

read_panel_genes <- function() {
  panel_path <- "results/tables/module_gene_panels_top20.tsv"
  if (file.exists(panel_path)) {
    panels <- readr::read_tsv(panel_path, show_col_types = FALSE)
    nm <- names(panels)

    panel_col <- intersect(c("panel", "Panel", "module", "Module", "score", "Score"), nm)[1]
    gene_col <- intersect(c("gene", "Gene", "gene_symbol", "symbol", "SYMBOL"), nm)[1]

    if (!is.na(panel_col) && !is.na(gene_col)) {
      p <- panels %>%
        mutate(panel_lower = stringr::str_to_lower(as.character(.data[[panel_col]])))

      prolif <- p %>%
        filter(str_detect(panel_lower, "prolif|panel_a|up")) %>%
        pull(.data[[gene_col]]) %>%
        unique()

      heploss <- p %>%
        filter(str_detect(panel_lower, "hep|panel_b|loss|down")) %>%
        pull(.data[[gene_col]]) %>%
        unique()

      if (length(prolif) >= 5 && length(heploss) >= 5) {
        return(list(ProlifHub = prolif, HepLoss = heploss))
      }
    }
  }

  # Fallback files from earlier pipeline versions.
  prolif_path <- "results/modules/prolifhub_genes_top20.txt"
  hep_path <- "results/modules/heploss_genes_top20.txt"

  if (file.exists(prolif_path) && file.exists(hep_path)) {
    return(list(
      ProlifHub = readr::read_lines(prolif_path) %>% discard(~ .x == ""),
      HepLoss = readr::read_lines(hep_path) %>% discard(~ .x == "")
    ))
  }

  stop(
    "Could not find module gene panels. Expected either:\n",
    "  results/tables/module_gene_panels_top20.tsv\n",
    "or:\n",
    "  results/modules/prolifhub_genes_top20.txt and results/modules/heploss_genes_top20.txt"
  )
}

# ============================================================
# Figure 5: HBV_INJURY coefficient plot
# ============================================================

make_fig3_hbv_injury_coefficients <- function() {
  summary_path <- "results/tables/gse121248_hbv_injury_topN_extended_manuscript_summary.tsv"
  fallback_path <- "results/tables/gse121248_cibersortx_adjusted_regression.tsv"

  if (file.exists(summary_path)) {
    coef <- readr::read_tsv(summary_path, show_col_types = FALSE)

    # Prefer TOP_2000 and EXTENDED_7792 if both exist.
    coef <- coef %>%
      mutate(
        injury_set = as.character(injury_set),
        analysis_role = as.character(analysis_role),
        model = as.character(model)
      )

    keep_sets <- c("HBV_INJURY_TOP_2000", "HBV_INJURY_EXTENDED_7792")
    if (any(coef$injury_set %in% keep_sets)) {
      coef <- coef %>% filter(injury_set %in% keep_sets)
    }

    coef_plot <- coef %>%
      transmute(
        injury_set,
        analysis_role,
        model,
        estimate = tumor_coefficient,
        conf.low = ci_low,
        conf.high = ci_high,
        p.value = p_value,
        percent_retained_vs_unadjusted
      )
  } else if (file.exists(fallback_path)) {
    coef <- readr::read_tsv(fallback_path, show_col_types = FALSE)
    coef_plot <- coef %>%
      mutate(injury_set = "HBV_INJURY", analysis_role = "primary") %>%
      rename(estimate = estimate, conf.low = conf.low, conf.high = conf.high, p.value = p.value) %>%
      mutate(percent_retained_vs_unadjusted = 100 * estimate / estimate[model == "unadjusted"]) %>%
      select(injury_set, analysis_role, model, estimate, conf.low, conf.high, p.value, percent_retained_vs_unadjusted)
  } else {
    stop("Missing HBV injury regression summary. Run script 09 first.")
  }

  # Normalize model labels across script versions.
  coef_plot <- coef_plot %>%
    mutate(
      model_short = case_when(
        model %in% c("unadjusted") ~ "Raw tumor effect",
        model %in% c("e2f_g2m_adjusted", "proliferation_adjusted") ~ "E2F/G2M adjusted",
        model %in% c("cibersortx_pc_adjusted", "proliferation_cibersortx_pc_adjusted") ~ "E2F/G2M + CIBERSORTx PCs",
        model %in% c("selected_fraction_adjusted", "proliferation_selected_fraction_adjusted") ~ "E2F/G2M + selected fractions",
        TRUE ~ model
      ),
      model_short = factor(
        model_short,
        levels = c(
          "Raw tumor effect",
          "E2F/G2M adjusted",
          "E2F/G2M + CIBERSORTx PCs",
          "E2F/G2M + selected fractions"
        )
      ),
      injury_set_label = case_when(
        injury_set == "HBV_INJURY_TOP_2000" ~ "HBV_INJURY_TOP_2000\ncompact primary",
        injury_set == "HBV_INJURY_EXTENDED_7792" ~ "HBV_INJURY_EXTENDED_7792\nextended sensitivity",
        TRUE ~ injury_set
      ),
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      p_text = paste0(sig_label, "\n", p_label(p.value))
    ) %>%
    filter(!is.na(model_short))

  readr::write_tsv(coef_plot, "results/tables/fig3_hbv_injury_coefficients_plot_data.tsv")

  fig3 <- ggplot(coef_plot, aes(x = estimate, y = model_short)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.55) +
    geom_point(size = 2.4) +
    geom_text(
      aes(label = p_text),
      nudge_y = 0.28,
      size = 3.0,
      lineheight = 0.95,
      check_overlap = TRUE
    ) +
    facet_wrap(~ injury_set_label, ncol = 1, scales = "free_y") +
    labs(
      x = "Tumor coefficient for HBV_INJURY score\n(tumor vs non-tumor/adjacent)",
      y = NULL,
      title = "Figure 5. HBV_INJURY tumor effect after proliferation and immune-composition adjustment",
      subtitle = "Points show tumor coefficients; horizontal bars show 95% confidence intervals"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggplot2::ggsave("results/figures/Fig5_HBV_INJURY_adjusted_coefficients.pdf", fig3, width = 7.2, height = 5.8)
  ggplot2::ggsave("results/figures/Fig5_HBV_INJURY_adjusted_coefficients.png", fig3, width = 7.2, height = 5.8, dpi = 400)
  ggplot2::ggsave("paper_package/figures/Fig5_HBV_INJURY_adjusted_coefficients.pdf", fig3, width = 7.2, height = 5.8)
  ggplot2::ggsave("paper_package/figures/Fig5_HBV_INJURY_adjusted_coefficients.png", fig3, width = 7.2, height = 5.8, dpi = 400)

  message("[Fig5] Wrote Fig5_HBV_INJURY_adjusted_coefficients to results/figures and paper_package/figures")
}

# ============================================================
# Figure 3: discovery-cohort module-score boxplots
# ============================================================

make_discovery_scores <- function(acc, panels) {
  obj_path <- file.path("data", "processed", paste0(acc, "_curated.rds"))
  if (!file.exists(obj_path)) stop("Missing curated object: ", obj_path)

  obj <- readr::read_rds(obj_path)
  expr <- obj$expr
  meta <- obj$metadata %>% mutate(sample_id = as.character(sample_id))

  if (!all(colnames(expr) %in% meta$sample_id)) {
    warning(acc, ": some expression columns are absent from metadata.")
  }

  tibble(
    sample_id = colnames(expr),
    ProlifHubScore = score_module_local(expr, panels$ProlifHub, min_genes = 5),
    HepLossScore = score_module_local(expr, panels$HepLoss, min_genes = 5)
  ) %>%
    mutate(HCCStateScore = ProlifHubScore - HepLossScore) %>%
    left_join(meta %>% select(sample_id, tissue), by = "sample_id") %>%
    mutate(
      dataset = acc,
      tissue = case_when(
        tissue %in% c("adjacent_or_nontumor", "adjacent", "normal", "non_tumor") ~ "non_tumor",
        tissue == "tumor" ~ "tumor",
        TRUE ~ as.character(tissue)
      ),
      tissue = factor(tissue, levels = c("non_tumor", "tumor"))
    ) %>%
    filter(tissue %in% c("non_tumor", "tumor"))
}

make_fig4_discovery_boxplots <- function() {
  panels <- read_panel_genes()

  scores <- bind_rows(
    make_discovery_scores("GSE121248", panels),
    make_discovery_scores("GSE41804", panels)
  )

  long <- scores %>%
    pivot_longer(
      cols = c(ProlifHubScore, HepLossScore, HCCStateScore),
      names_to = "score",
      values_to = "value"
    ) %>%
    mutate(
      score = factor(score, levels = c("ProlifHubScore", "HepLossScore", "HCCStateScore")),
      dataset = factor(dataset, levels = c("GSE121248", "GSE41804"))
    ) %>%
    filter(is.finite(value), !is.na(tissue))

  stats <- long %>%
    group_by(dataset, score) %>%
    summarise(
      n_tumor = sum(tissue == "tumor"),
      n_non_tumor = sum(tissue == "non_tumor"),
      delta = mean(value[tissue == "tumor"], na.rm = TRUE) - mean(value[tissue == "non_tumor"], na.rm = TRUE),
      cohen_d = safe_cohens_d(value[tissue == "tumor"], value[tissue == "non_tumor"]),
      p_value = tryCatch(stats::t.test(value ~ tissue)$p.value, error = function(e) NA_real_),
      y_pos = max(value, na.rm = TRUE) + 0.18 * diff(range(value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("p = ", ifelse(p_value < 1e-4, "<1e-4", signif(p_value, 2)),
                     "\nd = ", fmt_num(cohen_d, 2)),
      x_pos = 1.5,
      panel_id = paste0(LETTERS[as.integer(dataset)], as.integer(score))
    )

  readr::write_tsv(long, "results/tables/fig4_discovery_module_boxplot_data.tsv")
  readr::write_tsv(stats, "results/tables/fig4_discovery_module_boxplot_stats.tsv")

  # Panel labels: A/B for datasets, score panels via facet titles.
  dataset_labels <- c(
    GSE121248 = "A. GSE121248 (HBV-HCC discovery)",
    GSE41804 = "B. GSE41804 (HCV-HCC discovery)"
  )

  fig4 <- ggplot(long, aes(x = tissue, y = value)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.45) +
    geom_jitter(width = 0.12, height = 0, size = 0.85, alpha = 0.65) +
    geom_text(
      data = stats,
      aes(x = x_pos, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      lineheight = 0.95
    ) +
    facet_grid(
      rows = vars(dataset),
      cols = vars(score),
      scales = "free_y",
      labeller = labeller(dataset = dataset_labels)
    ) +
    labs(
      x = NULL,
      y = "Module score (mean gene-wise z-score)",
      title = "Figure 3. Discovery-cohort module score differences",
      subtitle = "Boxplots show ProlifHubScore, HepLossScore, and HCCStateScore in HBV- and HCV-associated HCC discovery cohorts"
    ) +
    scale_x_discrete(labels = c(non_tumor = "Non-tumor/\nadjacent", tumor = "Tumor")) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0),
      axis.text.x = element_text(size = 8.5),
      panel.grid.minor = element_blank()
    )

  ggplot2::ggsave("results/figures/Fig3_discovery_module_boxplots.pdf", fig4, width = 8.2, height = 5.8)
  ggplot2::ggsave("results/figures/Fig3_discovery_module_boxplots.png", fig4, width = 8.2, height = 5.8, dpi = 400)
  ggplot2::ggsave("paper_package/figures/Fig3_discovery_module_boxplots.pdf", fig4, width = 8.2, height = 5.8)
  ggplot2::ggsave("paper_package/figures/Fig3_discovery_module_boxplots.png", fig4, width = 8.2, height = 5.8, dpi = 400)

  message("[Fig3] Wrote Fig3_discovery_module_boxplots to results/figures and paper_package/figures")
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

make_fig3_hbv_injury_coefficients()
make_fig4_discovery_boxplots()

message("Done.")
