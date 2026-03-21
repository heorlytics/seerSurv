# seerSurv <img src="man/figures/logo.png" align="right" height="120" alt="seerSurv hex sticker" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/heorlytics/seerSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heorlytics/seerSurv/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/seerSurv)](https://CRAN.R-project.org/package=seerSurv)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**seerSurv** provides a transparent, reproducible framework for deriving monthly
cancer-recurrence **transition probabilities** used in state-transition health
economic models.

Given five-year aggregate survival proportions for **local/regional recurrence
(LR)** and **distant recurrence (DR)** populations—typically extracted from the
[SEER-17 registry](https://seer.cancer.gov/)—the package estimates:

| Symbol | Description |
|---|---|
| P\_LL | Monthly probability of remaining in LR |
| P\_DL | Monthly probability of progressing from LR to DR |
| P\_death given LR | Monthly probability of dying in LR (1 - P\_LL - P\_DL) |

Bundled data cover **11 cancer types** (Breast, Bladder, Colon & Rectum,
Esophageal, Kidney, Liver, Lung, Melanoma, Ovarian, Pancreatic, Stomach).

---

## Installation

```r
# Install the released CRAN version (once available)
install.packages("seerSurv")

# Install the development version from GitHub
# install.packages("pak")
pak::pak("heorlytics/seerSurv")
```

---

## Quick start

```r
library(seerSurv)

data(lifetable_seer)

result <- run_tumour_analysis(
  tumour       = "Melanoma",
  surv_vec_LR  = c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754),
  surv_vec_DR  = c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410),
  N_L          = 64775L,
  N_D          = 3121L,
  prop_male_LR = 0.580,
  mean_age_LR  = 61,
  prop_male_DR = 0.692,
  mean_age_DR  = 62,
  lifetable    = lifetable_seer,
  scenario     = "lifetime",   # "lifetime", "20y", or "5y"
  criterion    = "AIC"         # "AIC" or "BIC"
)

result
#> # A tibble: 1 x 7
#>     P_LL P_DL_approach2 P_Death_L_approach2 RMST_L_years RMST_D_years M_years M_months
#>    <dbl>          <dbl>               <dbl>        <dbl>        <dbl>   <dbl>    <dbl>
#> 1  0.999         0.0003               4e-04         17.2         3.15    14.1     169.
```

---

## Multi-tumour batch analysis

```r
library(dplyr)
library(tidyr)

data(tumour_data_seer)

results <- tumour_data_seer |>
  group_by(Tumor) |>
  summarise(
    out = list(
      run_tumour_analysis(
        tumour       = first(Tumor),
        surv_vec_LR  = c(S0[Recurrence == "LR"], S1[Recurrence == "LR"],
                         S2[Recurrence == "LR"], S3[Recurrence == "LR"],
                         S4[Recurrence == "LR"], S5[Recurrence == "LR"]),
        surv_vec_DR  = c(S0[Recurrence == "DR"], S1[Recurrence == "DR"],
                         S2[Recurrence == "DR"], S3[Recurrence == "DR"],
                         S4[Recurrence == "DR"], S5[Recurrence == "DR"]),
        N_L          = first(N_LR),
        N_D          = first(N_DR),
        prop_male_LR = prop_male[Recurrence == "LR"],
        mean_age_LR  = mean_age[Recurrence == "LR"],
        prop_male_DR = prop_male[Recurrence == "DR"],
        mean_age_DR  = mean_age[Recurrence == "DR"],
        lifetable    = lifetable_seer,
        scenario     = "lifetime",
        criterion    = "AIC"
      )
    ),
    .groups = "drop"
  ) |>
  unnest(out)
```

---

## Methodology

```
Input: 5-year aggregate survival (LR and DR) + demographics
     |
     v
 prep_ipd()           <- convert to pseudo-IPD (probability-weighted)
     |
     v
 fit_models()         <- 6 parametric distributions (flexsurv)
     |
     v
 extract_ic() +       <- select top-k by AIC / BIC
 compute_weights()    <- assign relative-likelihood weights
     |
     v
 blend_survival()     <- weighted-average survival curve
     |
     v
 make_background_surv() <- multiply by US-CDC life-table survival
     |
     v
 run_tumour_analysis()  <- optimise P(LL) and P(DL), return tibble
```

**P(LL)** is derived from the mean LR sojourn time (RMST differential).
**P(DL)** is estimated by minimising the relative sum-of-squares between the
observed LR survival curve and its convolution approximation.

Full details are provided in the companion vignette and the associated
publication.

---

## Bundled data

| Dataset | Description |
|---|---|
| `lifetable_seer` | US-CDC annual death rates by age and sex (ages 0-100) |
| `tumour_data_seer` | SEER-17 5-year survival proportions for 11 cancer types |

---

## Citation

If you use **seerSurv** in your research, please cite:

> Mansoori S, Pandey S, Rani R, Singh B, Kurt M (2025).
> "Deriving Cancer Progression Rates After Local/Regional Recurrence Using
> Aggregate Survival Data: An R Package for Health Economic Evaluations."
> *Value in Health* (submitted).

```r
citation("seerSurv")
```

---

## Contributing

Bug reports and pull requests are welcome on
[GitHub](https://github.com/heorlytics/seerSurv/issues).
Please read the [contribution guidelines](https://github.com/heorlytics/seerSurv/blob/master/CONTRIBUTING.md) before opening a
pull request.

---

## License

MIT © 2025 Sameer Mansoori, Shubhram Pandey, Rashi Rani, Barinder Singh,
Murat Kurt.  See [LICENSE](LICENSE) for details.
