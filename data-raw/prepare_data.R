## ============================================================================
## data-raw/prepare_data.R
##
## Reproducible script that reads the raw CSV files, applies light cleaning,
## and writes the binary .rda objects stored in data/.
##
## Run this script (once) with:
##   source("data-raw/prepare_data.R")
## or via devtools:
##   devtools::run_script("data-raw/prepare_data.R")
##
## This file is intentionally NOT part of the installed package.
## ============================================================================

library(tibble)
library(dplyr)

# ---- 1. lifetable_seer -------------------------------------------------------

lifetable_seer <- read.csv(
  file   = "data-raw/lifetable_seer.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Ensure Age is integer
lifetable_seer$Age <- as.integer(lifetable_seer$Age)

stopifnot(
  is.data.frame(lifetable_seer),
  all(c("Age", "Males", "Females") %in% names(lifetable_seer)),
  nrow(lifetable_seer) >= 100L,
  all(lifetable_seer$Males  >= 0 & lifetable_seer$Males  <= 1),
  all(lifetable_seer$Females >= 0 & lifetable_seer$Females <= 1)
)

usethis::use_data(lifetable_seer, overwrite = TRUE, compress = "xz")
message("lifetable_seer saved (", nrow(lifetable_seer), " rows)")


# ---- 2. tumour_data_seer -----------------------------------------------------
## This data frame is constructed directly from the SEER-17 extraction;
## source values are hard-coded here for full reproducibility.

tumour_data_seer <- tibble(
  Tumor = c(
    "Kidney",        "Kidney",
    "Bladder",       "Bladder",
    "Ovarian",       "Ovarian",
    "Breast",        "Breast",
    "Lung",          "Lung",
    "Liver",         "Liver",
    "Esophageal",    "Esophageal",
    "Colon & Rectum","Colon & Rectum",
    "Stomach",       "Stomach",
    "Pancreatic",    "Pancreatic",
    "Melanoma",      "Melanoma"
  ),
  Recurrence = rep(c("LR", "DR"), times = 11L),
  S0 = 1,
  S1 = c(
    0.9750, 0.4856,   0.8898, 0.3517,   0.9523, 0.7707,
    0.9960, 0.7313,   0.8156, 0.3589,   0.6352, 0.2112,
    0.7203, 0.3342,   0.9630, 0.6231,   0.8482, 0.3625,
    0.6709, 0.2461,   0.995947563, 0.549085165
  ),
  S2 = c(
    0.9559, 0.3427,   0.8103, 0.1943,   0.9103, 0.6272,
    0.9877, 0.5882,   0.7023, 0.2141,   0.4864, 0.0922,
    0.5474, 0.1549,   0.9334, 0.4237,   0.7267, 0.1774,
    0.4467, 0.0998,   0.988631992, 0.439560429
  ),
  S3 = c(
    0.9376, 0.2668,   0.7660, 0.1415,   0.8798, 0.5027,
    0.9779, 0.4780,   0.6285, 0.1544,   0.3979, 0.0572,
    0.4635, 0.0932,   0.9030, 0.2947,   0.6600, 0.1109,
    0.3488, 0.0605,   0.983051983, 0.394462897
  ),
  S4 = c(
    0.9212, 0.2245,   0.7359, 0.1189,   0.8476, 0.4132,
    0.9685, 0.3964,   0.5739, 0.1211,   0.3425, 0.0396,
    0.4147, 0.0725,   0.8754, 0.2197,   0.6213, 0.0841,
    0.2980, 0.0422,   0.978268696, 0.363009452
  ),
  S5 = c(
    0.9049, 0.1895,   0.7123, 0.1037,   0.8186, 0.3426,
    0.9596, 0.3339,   0.5317, 0.0978,   0.2997, 0.0281,
    0.3798, 0.0568,   0.8509, 0.1734,   0.5924, 0.0701,
    0.2720, 0.0342,   0.975359441, 0.340994938
  ),
  prop_male = c(
    0.649, 0.712,   0.790, 0.735,   0.000, 0.000,
    0.000, 0.000,   0.477, 0.552,   0.762, 0.723,
    0.804, 0.846,   0.565, 0.575,   0.612, 0.668,
    0.540, 0.573,   0.580, 0.692
  ),
  mean_age = c(
    61, 62,   65, 64,   59, 62,
    61, 61,   65, 64,   63, 63,
    63, 63,   61, 61,   63, 62,
    63, 63,   61, 62
  )
)

## Approximate SEER-17 sample sizes (informational)
n_by_tumor <- c(
  Kidney         = 12100, Bladder   = 7200,
  Ovarian        = 6800,  Breast    = 98000,
  Lung           = 18500, Liver     = 3100,
  Esophageal     = 2100, `Colon & Rectum` = 31000,
  Stomach        = 4900,  Pancreatic = 3400,
  Melanoma       = 67896
)
## Rough split: LR ~ 95 %, DR ~ 5 % (tumour-dependent; placeholder)
tumour_data_seer <- tumour_data_seer |>
  dplyr::left_join(
    tibble(
      Tumor = names(n_by_tumor),
      N_LR  = as.integer(n_by_tumor * 0.95),
      N_DR  = as.integer(n_by_tumor * 0.05)
    ),
    by = "Tumor"
  )

usethis::use_data(tumour_data_seer, overwrite = TRUE, compress = "xz")
message("tumour_data_seer saved (", nrow(tumour_data_seer), " rows)")
