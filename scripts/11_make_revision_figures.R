source("R/_shared.R")
suppressPackageStartupMessages({ library(patchwork); library(ggrepel) })
ensure_dirs()
# Fig 1: workflow placeholder generated as a clean ggplot diagram.
workflow <- tibble(x = c(1, 2.5, 4, 5.5, 7), y = 1,
                   label = c("A. Dataset selection\nGEO/TCGA registry\nInclusion/exclusion table",
                             "B. Discovery cohorts\nGSE121248 HBV-HCC\nGSE41804 HCV-HCC\nlimma + GSVA",
                             "C. Module construction\nMeta-analysis + I²/τ²\nProlifHub/HepLoss",
                             "D. External validation\nGSE14520 + added GEO cohorts\nAUC/effect sizes/null sets",
                             "E. Adjusted analyses\nTCGA clinical models\nESTIMATE/CIBERSORTx"))
p1 <- ggplot(workflow, aes(x, y, label = label)) +
  geom_label(size = 3.2, label.size = 0.3) +
  geom_segment(data = tibble(x = c(1.55, 3.05, 4.55, 6.05), xend = c(1.95, 3.45, 4.95, 6.45), y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "in"))) +
  theme_void() + xlim(0.3, 7.7) + ggtitle("Figure 1. Revised study design and analysis workflow")
dir.create("paper_package/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/Fig1_revised_workflow.pdf", p1, width = 11, height = 3.2)
ggsave("results/figures/Fig1_revised_workflow.png", p1, width = 11, height = 3.2, dpi = 400)
ggsave("paper_package/figures/Fig1_revised_workflow.pdf", p1, width = 11, height = 3.2)
ggsave("paper_package/figures/Fig1_revised_workflow.png", p1, width = 11, height = 3.2, dpi = 400)
# Fig 2: Hallmark dotplot.
hall <- readr::read_tsv("results/tables/discovery_hallmark_combined.tsv", show_col_types = FALSE) |>
  mutate(direction = ifelse(logFC > 0, "Activated in tumor", "Suppressed in tumor"), neglogFDR = -log10(FDR)) |>
  group_by(gene_set) |> mutate(best = min(FDR, na.rm = TRUE)) |> ungroup() |> arrange(best) |> slice_head(n = 40)
p2 <- ggplot(hall, aes(logFC, reorder(gene_set, logFC), size = neglogFDR, color = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(alpha = 0.85) + facet_wrap(~dataset) +
  labs(x = "Tumor vs non-tumor logFC", y = NULL, size = "-log10(FDR)", color = NULL) + theme_bw(base_size = 9)
ggsave("results/figures/Fig2_hallmark_dotplot.pdf", p2, width = 10, height = 8)
ggsave("results/figures/Fig2_hallmark_dotplot.png", p2, width = 10, height = 8, dpi = 400)
ggsave("paper_package/figures/Fig2_hallmark_dotplot.pdf", p2, width = 10, height = 8)
ggsave("paper_package/figures/Fig2_hallmark_dotplot.png", p2, width = 10, height = 8, dpi = 400)
# Fig 4: multicohort validation summary.
val <- readr::read_tsv("results/tables/geo_validation_summary.tsv", show_col_types = FALSE) |> filter(score == "HCCStateScore")
p5 <- ggplot(val, aes(x = delta_tumor_minus_nontumor, y = reorder(dataset, delta_tumor_minus_nontumor))) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(aes(size = AUC)) +
  geom_errorbarh(aes(xmin = delta_tumor_minus_nontumor - 1.96*abs(delta_tumor_minus_nontumor/cohens_d)/sqrt(n_tumor+n_non_tumor),
                     xmax = delta_tumor_minus_nontumor + 1.96*abs(delta_tumor_minus_nontumor/cohens_d)/sqrt(n_tumor+n_non_tumor)), height = 0.15) +
  labs(x = "HCCStateScore delta: tumor - non-tumor", y = NULL, size = "AUC") + theme_bw()
ggsave("results/figures/Fig4_multicohort_validation.pdf", p5, width = 7, height = 5)
ggsave("results/figures/Fig4_multicohort_validation.png", p5, width = 7, height = 5, dpi = 400)
dir.create("paper_package/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("paper_package/figures/Fig4_multicohort_validation.pdf", p5, width = 7, height = 5)
ggsave("paper_package/figures/Fig4_multicohort_validation.png", p5, width = 7, height = 5, dpi = 400)
# Fig robustness
if (file.exists("results/tables/module_size_robustness.tsv")) {
  rob <- readr::read_tsv("results/tables/module_size_robustness.tsv", show_col_types = FALSE)
  p <- ggplot(rob, aes(module_size, AUC, group = dataset)) + geom_line() + geom_point() + facet_wrap(~dataset) + theme_bw() + ylim(0.5, 1) +
    labs(x = "Top-N genes per module", y = "HCCStateScore AUC")
  ggsave("results/figures/SuppFig_module_size_robustness.pdf", p, width = 10, height = 7)
}
