# Estimate Cancer-Recurrence Transition Probabilities for One Tumour Type

The primary user-facing function of seerSurv. Given five-year aggregate
survival proportions for local/regional-recurrence (LR) and distant-
recurrence (DR) populations together with demographic inputs, this
function:

1.  Converts aggregate survival proportions to pseudo-IPD.

2.  Fits a panel of parametric models and selects the top-\\k\\ by
    information criterion (AIC or BIC).

3.  Blends the selected models using relative-likelihood weights.

4.  Adjusts for age- and sex-specific background mortality from US-CDC
    life-tables.

5.  Computes the mean LR sojourn time (RMST differential) and derives
    the monthly LR persistence probability \\P(\text{LR})\\.

6.  Estimates the monthly LR-to-DR transition probability
    \\P(\text{DR})\\ via a convolution optimisation.

7.  Returns all three monthly transition probabilities:
    \\P(\text{LR})\\, \\P(\text{DR})\\, and \\P(\text{Death} \|
    \text{LR})\\ = \\1 - P(\text{LR}) - P(\text{DR})\\.

## Usage

``` r
run_tumour_analysis(
  tumour,
  surv_vec_LR,
  surv_vec_DR,
  N_L,
  N_D,
  prop_male_LR,
  mean_age_LR,
  prop_male_DR,
  mean_age_DR,
  lifetable,
  scenario = c("lifetime", "20y", "5y"),
  criterion = c("AIC", "BIC"),
  dists = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma"),
  top_k = 3L,
  max_age = 100L
)
```

## Arguments

- tumour:

  Character scalar. Tumour label used for identification in multi-tumour
  analyses (e.g. `"Melanoma"`).

- surv_vec_LR:

  Numeric vector of survival probabilities for the LR population at
  years 0, 1, 2, 3, 4, 5 (must start at 1).

- surv_vec_DR:

  Numeric vector of survival probabilities for the DR population at
  years 0, 1, 2, 3, 4, 5 (must start at 1).

- N_L:

  Integer. Number of LR patients in the SEER cohort (informational;
  currently stored for provenance and future sample-size adjustments).

- N_D:

  Integer. Number of DR patients in the SEER cohort (informational).

- prop_male_LR:

  Numeric in \\\[0, 1\]\\. Proportion male in the LR population, used to
  construct a sex-mixed background mortality rate.

- mean_age_LR:

  Numeric. Mean (or median) age of the LR population.

- prop_male_DR:

  Numeric in \\\[0, 1\]\\. Proportion male in the DR population.

- mean_age_DR:

  Numeric. Mean (or median) age of the DR population.

- lifetable:

  A data frame with columns `Age`, `Males`, `Females` containing annual
  death rates. The bundled
  [`lifetable_seer`](https://heorlytics.github.io/seerSurv/reference/lifetable_seer.md)
  dataset is the canonical choice.

- scenario:

  Character scalar controlling the extrapolation horizon. One of:

  `"lifetime"`

  :   Lifetime horizon: \\100 - \min(\text{mean age LR}, \text{mean age
      DR})\\ years.

  `"20y"`

  :   Fixed 20-year horizon.

  `"5y"`

  :   Fixed 5-year horizon.

- criterion:

  Character scalar: `"AIC"` (default) or `"BIC"`. Governs both model
  selection and weight computation.

- dists:

  Character vector of distribution names passed to
  [`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md).
  Defaults to the six-distribution base set.

- top_k:

  Integer. Number of top-ranked models to include in the blend (default
  3).

- max_age:

  Integer. Maximum age assumed for the lifetime horizon calculation
  (default 100).

## Value

A one-row [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing:

- `P_LL`:

  Monthly probability of remaining in LR.

- `P_DL_approach2`:

  Monthly probability of transitioning from LR to DR (convolution
  optimisation, Approach 2).

- `P_Death_L_approach2`:

  Monthly probability of death from LR (\\1 - P\_{\text{LL}} -
  P\_{\text{DL}}\\).

- `RMST_L_years`:

  Net RMST for the LR population (years).

- `RMST_D_years`:

  Net RMST for the DR population (years).

- `M_years`:

  Mean LR sojourn time (RMST differential, years).

- `M_months`:

  Mean LR sojourn time (months).

## Details

### P(LR) derivation

The mean LR sojourn time \\M\\ (in months) is the area between the net
LR and net DR survival curves. \\P(\text{LR})\\ is the root of:
\$\$\frac{1 - p^{N+1}}{1 - p} = M,\$\$ where \\N\\ is the total number
of monthly time steps in the horizon.

### P(DR) derivation — Approach 2

The modelled LR survival \\R(t)\\ is approximated as the sum of patients
still in LR (\\L(t) = P\_{\text{LL}}^t\\) and those who have progressed
to DR: \$\$R(t) \approx L(t) + P\_{\text{DL}} \sum\_{i=0}^{t-1} L(i)\\
Q_D(t-i),\$\$ where \\Q_D(\cdot)\\ is the monthly DR survival.
\\P\_{\text{DL}}\\ is estimated by minimising the relative
sum-of-squares discrepancy between \\R(t)\\ and the right-hand side over
all \\t = 1, \ldots, N\\.

## See also

[`prep_ipd`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md),
[`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md),
[`compute_weights`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md),
[`blend_survival`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md),
[`make_background_surv`](https://heorlytics.github.io/seerSurv/reference/make_background_surv.md),
[`tumour_data_seer`](https://heorlytics.github.io/seerSurv/reference/tumour_data_seer.md)

## Examples

``` r
# \donttest{
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
  scenario     = "lifetime",
  criterion    = "AIC"
)
#> New names:
#> • `surv` -> `surv...1`
#> • `surv` -> `surv...2`
#> • `surv` -> `surv...3`
#> New names:
#> • `surv` -> `surv...1`
#> • `surv` -> `surv...2`
#> • `surv` -> `surv...3`
result
#> # A tibble: 1 × 7
#>    P_LL P_DL_approach2 P_Death_L_approach2 RMST_L_years RMST_D_years M_years
#>   <dbl>          <dbl>               <dbl>        <dbl>        <dbl>   <dbl>
#> 1 0.993        0.00659           0.0000467         21.2         9.18    12.0
#> # ℹ 1 more variable: M_months <dbl>
# }
```
