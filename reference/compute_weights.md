# Compute Relative-Likelihood Weights for Model Averaging

Selects the top-\\k\\ models by information criterion value and assigns
relative-likelihood weights using the Akaike / Schwarz weight formula:
\$\$w_i = \frac{\exp\\\bigl((\mathrm{IC}\_{\min} -
\mathrm{IC}\_i)/2\bigr)} {\sum\_{j=1}^{k}
\exp\\\bigl((\mathrm{IC}\_{\min} - \mathrm{IC}\_j)/2\bigr)}.\$\$

## Usage

``` r
compute_weights(ic_vals, top_k = 3L, criterion = "AIC")
```

## Arguments

- ic_vals:

  Named numeric vector of information criterion values, e.g. from
  [`extract_ic`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md).

- top_k:

  Integer. Number of top-ranked models to retain (default 3).

- criterion:

  Character scalar passed through for labelling only. Typically `"AIC"`
  or `"BIC"`.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
four columns:

- `model`:

  Distribution name (from names of `ic_vals`).

- `weight`:

  Relative-likelihood model weight (sums to 1).

- `ic`:

  Raw information criterion value.

- `criterion`:

  Label passed via the `criterion` argument.

## Examples

``` r
ic <- c(exp = 240.1, weibull = 238.7, gompertz = 241.5,
        llogis = 239.3, lnorm = 243.0, gamma = 238.9)
compute_weights(ic, top_k = 3, criterion = "AIC")
#> # A tibble: 3 × 4
#>   model   weight    ic criterion
#>   <chr>    <dbl> <dbl> <chr>    
#> 1 weibull  0.378  239. AIC      
#> 2 gamma    0.342  239. AIC      
#> 3 llogis   0.280  239. AIC      
```
