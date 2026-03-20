# Predict Survival from a Single `flexsurvreg` Model

A thin wrapper around
[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md)
that extracts the estimated survival probability at each requested time
point. Intended for internal use and for users who want individual model
curves before blending.

## Usage

``` r
predict_surv(model, times)
```

## Arguments

- model:

  A fitted
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  object.

- times:

  Numeric vector of time points at which to evaluate survival (same
  units as the data used to fit `model`).

## Value

Numeric vector of survival estimates, one per element of `times`.

## See also

[`blend_survival`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md)

## Examples

``` r
# \donttest{
t <- 0:5
s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd  <- prep_ipd(t, s)
mods <- fit_models(ipd[-1, ])
predict_surv(mods[["weibull"]], times = seq(0, 5, by = 0.5))
#>  [1] 1.0000000 0.9523380 0.9063679 0.8624061 0.8204459 0.7804346 0.7423041
#>  [8] 0.7059808 0.6713894 0.6384551 0.6071043
# }
```
