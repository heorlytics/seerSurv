# seerSurv: Cancer Progression Transition Probabilities from Aggregate Survival Data

The seerSurv package provides a transparent, reproducible framework for
estimating monthly transition probabilities used in state-transition
health economic models that include a local/regional recurrence (LR)
state.

Survival data for LR and distant-recurrence (DR) populations—typically
extracted from the SEER-17 registry—are extrapolated using a weighted
blend of parametric survival models. Age- and sex-specific background
mortality from US-CDC life-tables is layered on top to produce net
(all-cause) survival curves. The area differential between the LR and DR
curves yields the mean LR sojourn time, from which \\P(\text{LR})\\ is
derived. A convolution optimisation then recovers \\P(\text{DR})\\.

## Main workflow

1.  Prepare pseudo-individual patient data from aggregate survival
    proportions with
    [`prep_ipd`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md).

2.  Fit a panel of parametric models with
    [`fit_models`](https://heorlytics.github.io/seerSurv/reference/fit_models.md).

3.  Select and weight the top-\\k\\ models by AIC or BIC using
    [`extract_ic`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md)
    and
    [`compute_weights`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md).

4.  Blend the weighted survival curves with
    [`blend_survival`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md).

5.  Adjust for background mortality with
    [`make_background_surv`](https://heorlytics.github.io/seerSurv/reference/make_background_surv.md).

6.  Run the complete analysis pipeline for one or many tumour types with
    [`run_tumour_analysis`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md).

## Bundled data

- [`lifetable_seer`](https://heorlytics.github.io/seerSurv/reference/lifetable_seer.md):

  US-CDC annual death-rate life-table (ages 0–100) with separate columns
  for males and females, aligned to the SEER-17 study period.

- [`tumour_data_seer`](https://heorlytics.github.io/seerSurv/reference/tumour_data_seer.md):

  Reference 5-year aggregate survival proportions (years 0–5) for LR and
  DR populations across 11 cancer types derived from SEER-17, together
  with sample sizes, sex proportions, and mean ages.

## References

Mansoori S, Pandey S, Rani R, Singh B, Kurt M (2025). "Deriving Cancer
Progression Rates After Local/Regional Recurrence Using Aggregate
Survival Data: An R Package for Health Economic Evaluations." *Value in
Health* (submitted).

Surveillance, Epidemiology, and End Results (SEER) Program
(<https://seer.cancer.gov/>) SEER\*Stat Database.

## See also

Useful links:

- <https://github.com/heorlytics/seerSurv>

- <https://github.com/heorlytics/seerSurv/issues>

## Author

Sameer Mansoori, Shubhram Pandey, Rashi Rani, Barinder Singh, Murat Kurt

Maintainer: Shubhram Pandey <shubhram.pandey@heorlytics.com>
