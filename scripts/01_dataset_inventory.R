source("R/_shared.R")
ensure_dirs()
d <- read_datasets()
all_ds <- bind_rows(
  tibble(section = "discovery", accession = map_chr(d$discovery, "accession"), role = map_chr(d$discovery, "role"), etiology = map_chr(d$discovery, "etiology")),
  tibble(section = "hepatitis_axes", accession = map_chr(d$hepatitis_axes, "accession"), role = map_chr(d$hepatitis_axes, "role"), etiology = map_chr(d$hepatitis_axes, "etiology")),
  tibble(section = "validation", accession = map_chr(d$validation, "accession"), role = map_chr(d$validation, "role"), etiology = map_chr(d$validation, "etiology"))
) |> mutate(included = TRUE,
            inclusion_reason = case_when(
              section == "discovery" ~ "Paired or comparable tumor/non-tumor HCC cohort used for discovery.",
              section == "validation" ~ "Independent HCC cohort used for external validation.",
              TRUE ~ "Non-tumor hepatitis-stage cohort used to derive external biological axes."
            ),
            exclusion_reason = NA_character_)
readr::write_tsv(all_ds, "results/tables/dataset_inventory.tsv")
criteria <- tribble(
  ~criterion_type, ~criterion,
  "inclusion", "Human liver/HCC transcriptomic dataset with genome-wide expression profiling.",
  "inclusion", "Recoverable sample-level phenotype labels for tumor vs adjacent/non-tumor, or hepatitis-stage clinical strata for axis construction.",
  "inclusion", "At least 10 tumor and 10 non-tumor samples for tumor/non-tumor validation, unless used only for expression concordance or tumor-heavy sensitivity analysis.",
  "inclusion", "Processed matrices or count matrices available from GEO/TCGA with gene-level mapping possible.",
  "exclusion", "Cell-line-only, xenograft-only, or treatment-perturbation-only datasets without baseline liver tissue.",
  "exclusion", "Datasets lacking usable sample-level phenotype labels after manual metadata inspection.",
  "exclusion", "Datasets with too few samples for group-level inference.",
  "exclusion", "Duplicated cohorts or subsets already represented by a larger selected dataset."
)
readr::write_tsv(criteria, "results/tables/inclusion_exclusion_criteria.tsv")
