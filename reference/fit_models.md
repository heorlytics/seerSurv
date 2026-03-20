# Fit a Panel of Parametric Survival Models to Pseudo-IPD

Fits each distribution named in `dists` to a pseudo-individual patient
data (pseudo-IPD) data frame using
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
with observation weights. The resulting named list is designed to be
passed directly to
[`extract_ic`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md)
and
[`blend_survival`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md).

## Usage

``` r
fit_models(
  ipd,
  dists = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")
)
```

## Arguments

- ipd:

  A tibble / data frame with at least three columns:

  `time`

  :   Event or censoring time (numeric, \> 0).

  `event`

  :   Event indicator (1 = event, 0 = censored).

  `weight`

  :   Observation weight (non-negative numeric).

  Typically the output of
  [`prep_ipd`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md)
  with the first row (time = 0) removed.

- dists:

  Character vector of distribution names accepted by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  Defaults to the six-distribution set used throughout seerSurv:
  `c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")`.

## Value

A named list of `flexsurvreg` objects, one per distribution. Names match
the corresponding element of `dists`.

## Details

The `ipd` data frame should have the first row (time = 0) removed before
calling this function because `flexsurvreg` requires strictly positive
event times. A common idiom is:

      ipd  <- prep_ipd(time_vec, surv_vec)
      mods <- fit_models(ipd[-1, ])

## See also

[`prep_ipd`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md),
[`extract_ic`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md),
[`blend_survival`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md)

## Examples

``` r
# \donttest{
t <- 0:5
s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)   # melanoma DR (SEER-17)
ipd   <- prep_ipd(t, s)
mods  <- fit_models(ipd[-1, ])
names(mods)
#> [1] "exp"      "weibull"  "gompertz" "llogis"   "lnorm"    "gamma"   
# }
```
