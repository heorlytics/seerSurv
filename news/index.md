# Changelog

## seerSurv 0.1.0 (2026-03-21)

### Initial release

#### New features

- [`run_tumour_analysis()`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md):
  Primary wrapper function. Estimates monthly LR transition
  probabilities from aggregate survival data. Supports three time
  horizons (`"lifetime"`, `"20y"`, `"5y"`) and two model-selection
  criteria (`"AIC"`, `"BIC"`).

- [`prep_ipd()`](https://heorlytics.github.io/seerSurv/reference/prep_ipd.md):
  Converts aggregate survival proportions to probability-weighted
  pseudo-individual patient data.

- [`fit_models()`](https://heorlytics.github.io/seerSurv/reference/fit_models.md):
  Fits a panel of six parametric distributions (exponential, Weibull,
  Gompertz, log-logistic, log-normal, gamma) to pseudo-IPD via
  `flexsurv`.

- [`extract_ic()`](https://heorlytics.github.io/seerSurv/reference/extract_ic.md):
  Extracts AIC or BIC from a named list of `flexsurvreg` models.

- [`compute_weights()`](https://heorlytics.github.io/seerSurv/reference/compute_weights.md):
  Selects the top-k models and assigns relative-likelihood
  model-averaging weights.

- [`blend_survival()`](https://heorlytics.github.io/seerSurv/reference/blend_survival.md):
  Produces a weighted-average blended survival curve.

- [`predict_surv()`](https://heorlytics.github.io/seerSurv/reference/predict_surv.md):
  Wrapper to predict survival at arbitrary time points from a single
  `flexsurvreg` model.

- [`make_background_surv()`](https://heorlytics.github.io/seerSurv/reference/make_background_surv.md):
  Constructs background (all-cause) survival from age- and sex-specific
  US-CDC life-table rates.

- [`rmst()`](https://heorlytics.github.io/seerSurv/reference/rmst.md):
  Restricted mean survival time (trapezoidal rule).

- [`plot_survival_curves()`](https://heorlytics.github.io/seerSurv/reference/plot_survival_curves.md):
  `ggplot2`-based visualisation of up to four survival curves (blended
  LR, blended DR, net LR, net DR).

#### Bundled data

- `lifetable_seer`: US-CDC annual death rates by single year of age and
  sex (ages 0–100), aligned to the SEER-17 study period.

- `tumour_data_seer`: SEER-17 five-year aggregate survival proportions
  for 11 cancer types (Breast, Bladder, Colon & Rectum, Esophageal,
  Kidney, Liver, Lung, Melanoma, Ovarian, Pancreatic, Stomach).
