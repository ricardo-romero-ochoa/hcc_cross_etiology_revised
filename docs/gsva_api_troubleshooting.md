# GSVA API troubleshooting

Recent Bioconductor versions of `GSVA` no longer support calling `GSVA::gsva()` directly on a matrix. The old call:

```r
GSVA::gsva(expr, sets, method = "gsva", kcdf = "Gaussian")
```

can fail with:

```text
unable to find an inherited method for function 'gsva' for signature 'param = "matrix"'
```

The revision repo now uses `run_gsva_scores()` from `R/_shared.R`, which first tries the new `gsvaParam()` API and then falls back to the old matrix API only when appropriate.

The patched `scripts/03_discovery_hallmark_gsva.R` now calls:

```r
gs <- run_gsva_scores(expr, sets, method = "gsva", kcdf = "Gaussian",
                      min_size = cfg$gsva$min_size %||% 10,
                      max_size = cfg$gsva$max_size %||% 500,
                      verbose = FALSE)
```
