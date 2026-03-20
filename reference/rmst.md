# Compute the Restricted Mean Survival Time (Trapezoidal Rule)

Calculates the area under a survival curve up to the last observed time
point using the composite trapezoidal rule. The result equals the
restricted mean survival time (RMST) when the input curve has been
evaluated on a sufficiently fine grid.

## Usage

``` r
rmst(time, surv)
```

## Arguments

- time:

  Numeric vector of time points (strictly increasing, same units as the
  survival data).

- surv:

  Numeric vector of survival probabilities evaluated at each element of
  `time`. Must be the same length as `time`.

## Value

A single numeric value: the area under the curve (RMST) in the same
units as `time`.

## Details

The trapezoidal approximation is: \$\$\text{RMST} \approx
\sum\_{i=1}^{n-1} (t\_{i+1} - t_i) \cdot \frac{S(t_i) +
S(t\_{i+1})}{2}.\$\$ For uniform grids this simplifies to: \$\$\Delta t
\cdot \Bigl\[\tfrac{1}{2}S(t_0) + S(t_1) + \cdots + S(t\_{n-2}) +
\tfrac{1}{2}S(t\_{n-1})\Bigr\].\$\$

## See also

[`run_tumour_analysis`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md)

## Examples

``` r
# RMST from a synthetic exponential curve
t <- seq(0, 10, by = 0.1)
s <- exp(-0.1 * t)
rmst(t, s)   # close to 1/0.1 = 10 as horizon -> infinity
#> [1] 6.321258
```
