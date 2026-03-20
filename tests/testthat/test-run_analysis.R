## ============================================================================
## tests/testthat/test-run_analysis.R
##
## Integration tests for run_tumour_analysis().  These are slow; skip on CRAN.
## ============================================================================

skip_on_cran()
skip_if_not_installed("flexsurv")

data(lifetable_seer, package = "seerTP", envir = environment())

s_lr <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
s_dr <- c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410)

# ---- helper: run melanoma with default args --------------------------------- #
melanoma_result <- function(scenario = "lifetime", criterion = "AIC") {
  run_tumour_analysis(
    tumour       = "Melanoma",
    surv_vec_LR  = s_lr,
    surv_vec_DR  = s_dr,
    N_L          = 64775L,
    N_D          = 3121L,
    prop_male_LR = 0.580,
    mean_age_LR  = 61,
    prop_male_DR = 0.692,
    mean_age_DR  = 62,
    lifetable    = lifetable_seer,
    scenario     = scenario,
    criterion    = criterion
  )
}

# ---- structure -------------------------------------------------------------- #
test_that("run_tumour_analysis returns a one-row tibble", {
  out <- melanoma_result()
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
})

test_that("run_tumour_analysis returns all expected columns", {
  out <- melanoma_result()
  expected_cols <- c(
    "P_LL", "P_DL_approach2", "P_Death_L_approach2",
    "RMST_L_years", "RMST_D_years", "M_years", "M_months"
  )
  expect_named(out, expected_cols)
})

# ---- probability constraints ------------------------------------------------ #
test_that("P_LL lies in (0, 1)", {
  out <- melanoma_result()
  expect_gt(out$P_LL, 0)
  expect_lt(out$P_LL, 1)
})

test_that("P_DL_approach2 is non-negative", {
  out <- melanoma_result()
  expect_gte(out$P_DL_approach2, 0)
})

test_that("P_Death_L_approach2 is non-negative", {
  out <- melanoma_result()
  expect_gte(out$P_Death_L_approach2, 0)
})

test_that("monthly probabilities sum to <= 1", {
  out <- melanoma_result()
  total <- out$P_LL + out$P_DL_approach2 + out$P_Death_L_approach2
  expect_lte(total, 1 + 1e-6)
})

# ---- RMST ordering ---------------------------------------------------------- #
test_that("LR RMST > DR RMST (LR patients survive longer)", {
  out <- melanoma_result()
  expect_gt(out$RMST_L_years, out$RMST_D_years)
})

test_that("sojourn time M is positive", {
  out <- melanoma_result()
  expect_gt(out$M_years, 0)
  expect_gt(out$M_months, 0)
})

# ---- scenario sensitivity --------------------------------------------------- #
test_that("5-year horizon gives smaller RMST than lifetime", {
  r_life <- melanoma_result("lifetime")
  r_5y   <- melanoma_result("5y")
  expect_lt(r_5y$RMST_L_years, r_life$RMST_L_years)
})

test_that("BIC and AIC results differ (different weights)", {
  r_aic <- melanoma_result("lifetime", "AIC")
  r_bic <- melanoma_result("lifetime", "BIC")
  # Values should not be identical
  expect_false(isTRUE(all.equal(r_aic$P_LL, r_bic$P_LL, tolerance = 1e-10)))
})

# ---- input validation ------------------------------------------------------- #
test_that("run_tumour_analysis errors on invalid prop_male", {
  expect_error(
    run_tumour_analysis(
      tumour = "X", surv_vec_LR = s_lr, surv_vec_DR = s_dr,
      N_L = 1, N_D = 1,
      prop_male_LR = 1.5, mean_age_LR = 61,    # invalid
      prop_male_DR = 0.5,  mean_age_DR = 62,
      lifetable = lifetable_seer, scenario = "5y", criterion = "AIC"
    ),
    "\\[0, 1\\]"
  )
})

test_that("run_tumour_analysis errors on missing lifetable columns", {
  bad_lt <- lifetable_seer[, c("Age", "Males")]   # missing Females
  expect_error(
    run_tumour_analysis(
      tumour = "X", surv_vec_LR = s_lr, surv_vec_DR = s_dr,
      N_L = 1, N_D = 1,
      prop_male_LR = 0.5, mean_age_LR = 61,
      prop_male_DR = 0.5, mean_age_DR = 62,
      lifetable = bad_lt, scenario = "5y", criterion = "AIC"
    ),
    "Females"
  )
})
