# Blend Parametric Survival Curves Using Model-Averaging Weights

Produces a single, weighted-average survival curve by combining
predicted survival probabilities from the top-\\k\\ parametric models
selected by
[`compute_weights`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md).
The blend is a convex combination: \$\$\hat{S}(t) = \sum\_{i=1}^{k} w_i
\\ S_i(t),\$\$ where \\S_i(t)\\ is the predicted survival from the
\\i\\-th model and \\w_i\\ are the relative-likelihood weights.

## Usage

``` r
blend_survival(models, weights_tbl, times)
```

## Arguments

- models:

  A named list of fitted
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  objects, typically from
  [`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md).

- weights_tbl:

  A tibble from
  [`compute_weights`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md)
  with columns `model` (matching names in `models`) and `weight`.

- times:

  Numeric vector of time points at which to evaluate the blended
  survival curve.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
two columns:

- `time`:

  Time points (same as `times`).

- `surv`:

  Blended survival probability.

## Examples

``` r
# \donttest{
t <- 0:5
s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd  <- prep_ipd(t, s)
mods <- fit_models(ipd[-1, ])
ic   <- extract_ic(mods, "AIC")
wts  <- compute_weights(ic, top_k = 3, criterion = "AIC")
grid <- seq(0, 10, by = 0.25)
S_blend <- blend_survival(mods, wts, grid)
#> New names:
#> • `surv` -> `surv...1`
#> • `surv` -> `surv...2`
#> • `surv` -> `surv...3`
head(S_blend)
#> # A tibble: 6 × 2
#>    time  surv
#>   <dbl> <dbl>
#> 1  0    1    
#> 2  0.25 0.976
#> 3  0.5  0.951
#> 4  0.75 0.925
#> 5  1    0.901
#> 6  1.25 0.876
# }
```
