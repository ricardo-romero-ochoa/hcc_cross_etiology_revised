options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("dplyr", "tidyr", "tibble", "stringr", "readr", "purrr", "yaml", "ggplot2", "pROC",
          "metafor", "survival", "survminer", "ggrepel", "patchwork", "data.table", "matrixStats", "broom")
bioc <- c("GEOquery", "Biobase", "limma", "GSVA", "msigdbr", "SummarizedExperiment", "TCGAbiolinks", "estimate", "edgeR")
missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran, dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
