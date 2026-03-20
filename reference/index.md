# Package index

## Main analysis wrapper

End-to-end function for a single tumour type

- [`run_tumour_analysis()`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md)
  : Estimate Cancer-Recurrence Transition Probabilities for One Tumour
  Type

## Step-by-step building blocks

Individual functions that run_tumour_analysis() calls internally,
exposed for custom workflows.

- [`prep_ipd()`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md)
  : Prepare Pseudo-Individual Patient Data from Aggregate Survival
  Proportions

- [`fit_models()`](https://heorlytics.github.io/seerSurv/reference/fit_models.md)
  : Fit a Panel of Parametric Survival Models to Pseudo-IPD

- [`extract_ic()`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md)
  : Extract Information Criteria from a List of Parametric Survival
  Models

- [`compute_weights()`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md)
  : Compute Relative-Likelihood Weights for Model Averaging

- [`blend_survival()`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md)
  : Blend Parametric Survival Curves Using Model-Averaging Weights

- [`predict_surv()`](https://heorlytics.github.io/seerSurv/reference/predict_surv.md)
  :

  Predict Survival from a Single `flexsurvreg` Model

- [`make_background_surv()`](https://heorlytics.github.io/seerSurv/reference/make_background_surv.md)
  : Construct Background (All-Cause) Survival from a US-CDC Life-Table

- [`rmst()`](https://heorlytics.github.io/seerSurv/reference/rmst.md) :
  Compute the Restricted Mean Survival Time (Trapezoidal Rule)

## Visualisation

- [`plot_survival_curves()`](https://heorlytics.github.io/seerSurv/reference/plot_survival_curves.md)
  : Plot Blended Survival Curves for LR and DR Populations

## Bundled datasets

- [`lifetable_seer`](https://heorlytics.github.io/seerSurv/reference/lifetable_seer.md)
  : US-CDC Life-Table for SEER-17 Background Mortality Adjustment
- [`tumour_data_seer`](https://heorlytics.github.io/seerSurv/reference/tumour_data_seer.md)
  : SEER-17 Five-Year Aggregate Survival Data for 11 Cancer Types

## Package

- [`seerSurv`](https://heorlytics.github.io/seerSurv/reference/seerSurv-package.md)
  [`seerSurv-package`](https://heorlytics.github.io/seerSurv/reference/seerSurv-package.md)
  : seerSurv: Cancer Progression Transition Probabilities from Aggregate
  Survival Data
