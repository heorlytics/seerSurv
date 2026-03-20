# Prepare Pseudo-Individual Patient Data from Aggregate Survival Proportions

Converts a vector of aggregate survival proportions (typically read from
a Kaplan–Meier table or a registry publication) into a pseudo-individual
patient data (pseudo-IPD) tibble suitable for passing to
[`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md).
Each row represents one "pseudo-patient"; its weight is the probability
mass that the Kaplan–Meier curve drops at that time point. The last
observed subject is treated as censored.

## Usage

``` r
prep_ipd(time, surv)
```

## Arguments

- time:

  Numeric vector of observed time points (in years). Must be
  non-negative, strictly increasing, and the same length as `surv`.

- surv:

  Numeric vector of survival probabilities at each time point. Must
  start at 1 (or very close to it), be non-increasing, and lie in \\\[0,
  1\]\\.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
three columns:

- `time`:

  Time point (same units as input).

- `weight`:

  Probability mass at that time point (non-negative, sums to 1).

- `event`:

  Event indicator: `1` for an observed event, `0` for the last
  (censored) observation.

Rows with zero weight are silently dropped.

## Details

The weight for time \\t_i\\ is computed as \$\$w_i = \| S(t_i) -
S(t\_{i+1}) \|,\$\$ with the final weight adjusted so that all weights
sum to 1. This ensures that the empirical distribution implied by the
pseudo-IPD exactly reproduces the input survival curve.

## Examples

``` r
# Five-year melanoma LR survival (SEER-17)
t <- 0:5
s <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
ipd <- prep_ipd(t, s)
ipd
#> # A tibble: 6 × 3
#>    time  weight event
#>   <int>   <dbl> <int>
#> 1     0 0.00410     1
#> 2     1 0.00730     1
#> 3     2 0.00550     1
#> 4     3 0.00480     1
#> 5     4 0.00290     1
#> 6     5 0.975       0
```
