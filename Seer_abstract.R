library(dplyr)
library(tibble)
library(survival)
library(flexsurv)
library(purrr)
library(dplyr)
library(ggplot2)
library(readxl)

options(scipen = 99999999)

source("helper funtions/helper_functions.R")

lifetable <- read.csv("data/lifetable_seer.csv")
tumour_df <- read_xlsx("data/tumour_data.xlsx")

dists <- c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")
horizon  = "lifetime" #lifetime/1-100
criterion = "AIC" # AIC/BIC

melanoma_df <- tumour_df %>% filter(Tumor == "Melanoma")

results_melanoma <- melanoma_df %>%
  group_by(Tumor) %>%
  summarise(
    out = list(
      run_tumour_analysis(
        tumour = first(Tumor),
        surv_vec_LR = c(S0[Recurrence == "LR"], S1[Recurrence == "LR"],
                        S2[Recurrence == "LR"], S3[Recurrence == "LR"],
                        S4[Recurrence == "LR"], S5[Recurrence == "LR"]),
        surv_vec_DR = c(S0[Recurrence == "DR"], S1[Recurrence == "DR"],
                        S2[Recurrence == "DR"], S3[Recurrence == "DR"],
                        S4[Recurrence == "DR"], S5[Recurrence == "DR"]),
        N_L = N[Recurrence == "LR"],
        N_D = N[Recurrence == "DR"],
        prop_male_LR = prop_male[Recurrence == "LR"],
        mean_age_LR  = mean_age[Recurrence == "LR"],
        prop_male_DR = prop_male[Recurrence == "DR"],
        mean_age_DR  = mean_age[Recurrence == "DR"],
        lifetable    = lifetable,
        horizon      = horizon,
        criterion    = criterion
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest(out)

print(as.data.frame(results_melanoma$out$results), digits = 6)
results_melanoma$out$plots


# results_all <- tumour_df %>%
#   group_by(Tumor) %>%
#   summarise(
#     out = list(
#       run_tumour_analysis(
#         tumour = first(Tumor),
#         surv_vec_LR = c(S0[Recurrence == "LR"], S1[Recurrence == "LR"],
#                         S2[Recurrence == "LR"], S3[Recurrence == "LR"],
#                         S4[Recurrence == "LR"], S5[Recurrence == "LR"]),
#         surv_vec_DR = c(S0[Recurrence == "DR"], S1[Recurrence == "DR"],
#                         S2[Recurrence == "DR"], S3[Recurrence == "DR"],
#                         S4[Recurrence == "DR"], S5[Recurrence == "DR"]),
#         N_L = N[Recurrence == "LR"],
#         N_D = N[Recurrence == "DR"],
#         prop_male_LR = prop_male[Recurrence == "LR"],
#         mean_age_LR  = mean_age[Recurrence == "LR"],
#         prop_male_DR = prop_male[Recurrence == "DR"],
#         mean_age_DR  = mean_age[Recurrence == "DR"],
#         lifetable = lifetable,
#         horizon = horizon,
#         criterion = criterion
#       )
#     ),
#     .groups = "drop"
#   ) %>%
#   tidyr::unnest(out)
#
# print(as.data.frame(results_all$out$results), digits = 6)


# results_raw <- tumour_df %>%
#   group_by(Tumor) %>%
#   summarise(
#     out = list(
#       run_tumour_analysis(
#         tumour       = first(Tumor),
#         surv_vec_LR  = c(S0[Recurrence=="LR"], S1[Recurrence=="LR"],
#                          S2[Recurrence=="LR"], S3[Recurrence=="LR"],
#                          S4[Recurrence=="LR"], S5[Recurrence=="LR"]),
#         surv_vec_DR  = c(S0[Recurrence=="DR"], S1[Recurrence=="DR"],
#                          S2[Recurrence=="DR"], S3[Recurrence=="DR"],
#                          S4[Recurrence=="DR"], S5[Recurrence=="DR"]),
#         N_L          = N[Recurrence=="LR"],
#         N_D          = N[Recurrence=="DR"],
#         prop_male_LR = prop_male[Recurrence=="LR"],
#         mean_age_LR  = mean_age[Recurrence=="LR"],
#         prop_male_DR = prop_male[Recurrence=="DR"],
#         mean_age_DR  = mean_age[Recurrence=="DR"],
#         lifetable    = lifetable,
#         horizon      = horizon,
#         criterion    = criterion
#       )
#     ),
#     .groups = "drop"
#   )
#
# results_table <- dplyr::bind_rows(
#   lapply(results_raw$out, function(x) x$results)
# ) %>%
#   dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
#   dplyr::select(
#     Tumour        = tumour,
#     Horizon       = horizon_label,
#     Criterion     = criterion,
#     `P(L|L)`      = P_LL,
#     `P(D|L)`      = P_DL_approach2,
#     `P(Death|L)`  = P_Death_L_approach2,
#     `M (months)`  = M_months,
#     `RMST LR (y)` = RMST_L_years,
#     `RMST DR (y)` = RMST_D_years
#   )
#
# print(results_table, n = Inf)
# View(results_table)
