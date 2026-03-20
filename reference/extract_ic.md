# Extract Information Criteria from a List of Parametric Survival Models

Extracts either the Akaike Information Criterion (AIC) or the Bayesian
Information Criterion (BIC) from each `flexsurvreg` object in a named
list, returning a named numeric vector for use with
[`compute_weights`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md).

## Usage

``` r
extract_ic(models, criterion = c("AIC", "BIC"))
```

## Arguments

- models:

  A named list of fitted
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  objects, typically the output of
  [`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md).

- criterion:

  Character scalar: `"AIC"` (default) or `"BIC"`.

## Value

A named numeric vector of information criterion values, one per model.

## Examples

``` r
# \donttest{
t <- 0:5
s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd   <- prep_ipd(t, s)
mods  <- fit_models(ipd[-1, ])
extract_ic(mods, "AIC")
#>      exp  weibull gompertz   llogis    lnorm    gamma 
#> 3.375273 5.375252 5.356812 5.357385 5.331519 5.374326 
extract_ic(mods, "BIC")
#>       exp   weibull  gompertz    llogis     lnorm     gamma 
#> 0.7756162 0.1759382 0.1574987 0.1580710 0.1322055 0.1750123 
# }
```
