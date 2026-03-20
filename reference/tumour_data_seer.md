# SEER-17 Five-Year Aggregate Survival Data for 11 Cancer Types

Reference aggregate survival proportions at years 0–5 post-recurrence
for both local/regional-recurrence (LR) and distant-recurrence (DR)
populations across 11 cancer types, extracted from the Surveillance,
Epidemiology, and End Results SEER-17 registry (cases diagnosed before
2021). Also includes sample sizes, proportion male, and mean age for
each subgroup.

## Usage

``` r
tumour_data_seer
```

## Format

A data frame with 22 rows (2 recurrence types × 11 tumours) and 12
columns:

- `Tumor`:

  Cancer type label.

- `Recurrence`:

  Recurrence type: `"LR"` or `"DR"`.

- `S0`:

  Survival at year 0 (always 1).

- `S1`:

  Survival at year 1.

- `S2`:

  Survival at year 2.

- `S3`:

  Survival at year 3.

- `S4`:

  Survival at year 4.

- `S5`:

  Survival at year 5.

- `prop_male`:

  Proportion of patients who are male.

- `mean_age`:

  Mean (or median) age at recurrence.

- `N_LR`:

  Number of LR patients in the SEER cohort.

- `N_DR`:

  Number of DR patients in the SEER cohort.

## Source

Surveillance, Epidemiology, and End Results (SEER) Program
(<https://seer.cancer.gov/>) SEER\*Stat Database: SEER 17 Regs Research
Data, Nov 2022 Sub (2000–2020).

## Details

The 11 tumour types are: Breast, Bladder, Colon & Rectum, Esophageal,
Kidney, Liver, Lung, Melanoma, Ovarian, Pancreatic, and Stomach.

## See also

[`run_tumour_analysis`](https://heorlytics.github.io/seerSurv/reference/run_tumour_analysis.md)

## Examples

``` r
data(tumour_data_seer)
dplyr::filter(tumour_data_seer, Tumor == "Melanoma")
#> # A tibble: 2 × 12
#>   Tumor  Recurrence    S0    S1    S2    S3    S4    S5 prop_male mean_age  N_LR
#>   <chr>  <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>    <dbl> <int>
#> 1 Melan… LR             1 0.996 0.989 0.983 0.978 0.975     0.58        61 64501
#> 2 Melan… DR             1 0.549 0.440 0.394 0.363 0.341     0.692       62 64501
#> # ℹ 1 more variable: N_DR <int>
```
