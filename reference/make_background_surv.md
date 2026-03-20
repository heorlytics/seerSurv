# Construct Background (All-Cause) Survival from a US-CDC Life-Table

Uses age- and sex-specific annual death rates from a life-table to build
a cumulative background survival curve, then interpolates it onto a
user- supplied time grid. The resulting curve is multiplied against a
relative (disease-specific) survival curve inside
[`run_tumour_analysis`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md)
to obtain net (all-cause) survival.

## Usage

``` r
make_background_surv(lifetable, mean_age, rate_col, time_grid_years)
```

## Arguments

- lifetable:

  A data frame containing the life-table. Must include an `Age` column
  (integer ages 0, 1, 2, …) and at least one numeric column whose name
  is supplied via `rate_col`. The bundled
  [`lifetable_seer`](https://heorlytics.github.io/seerSurv/reference/lifetable_seer.md)
  dataset is the canonical input.

- mean_age:

  Numeric scalar. The mean (or median) age of the study cohort. The
  life-table is sub-set to rows where `Age >= floor(mean_age)`.

- rate_col:

  Character scalar. Name of the column in `lifetable` that contains the
  annual death rates to use (e.g. `"btrate_LR"` or `"btrate_DR"` after
  sex-mixing has been applied).

- time_grid_years:

  Numeric vector of time points (in years) onto which the cumulative
  survival is interpolated.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
two columns:

- `time`:

  Time in years (same as `time_grid_years`).

- `surv`:

  Background survival probability at each time point.

## Details

The cumulative background survival at integer year \\t\\ (starting from
`floor(mean_age)`) is: \$\$B(t) = \prod\_{j=0}^{t-1} \bigl(1 -
q_j\bigr),\$\$ where \\q_j\\ is the annual death rate at age \\j +
\lfloor\text{mean age}\rfloor\\. Non-integer time points are linearly
interpolated; times beyond the last life-table row use flat
extrapolation (`rule = 2` in
[`approx`](https://rdrr.io/r/stats/approxfun.html)).

## See also

[`lifetable_seer`](https://heorlytics.github.io/seerSurv/reference/lifetable_seer.md),
[`run_tumour_analysis`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md)

## Examples

``` r
data(lifetable_seer)

# Add a sex-mixed rate column (60% male cohort)
lt <- lifetable_seer
lt$btrate_mix <- lt$Males * 0.60 + lt$Females * 0.40

grid <- seq(0, 39, by = 1)   # 39-year horizon from age 61
B <- make_background_surv(lt, mean_age = 61, rate_col = "btrate_mix",
                          time_grid_years = grid)
head(B)
#> # A tibble: 6 × 2
#>    time  surv
#>   <dbl> <dbl>
#> 1     0 1    
#> 2     1 0.990
#> 3     2 0.978
#> 4     3 0.967
#> 5     4 0.954
#> 6     5 0.941
```
