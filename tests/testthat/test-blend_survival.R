## ============================================================================
## tests/testthat/test-blend_survival.R
## ============================================================================

skip_if_not_installed("flexsurv")

# Shared fixtures ------------------------------------------------------------- #
t_dr  <- 0:5
s_dr  <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd_D <- prep_ipd(t_dr, s_dr)
mods  <- fit_models(ipd_D[-1, ])
ic    <- extract_ic(mods, "AIC")
wts   <- compute_weights(ic, top_k = 3, criterion = "AIC")
grid  <- seq(0, 10, by = 1)
# ----------------------------------------------------------------------------- #

test_that("blend_survival returns a tibble with columns time and surv", {
  out <- blend_survival(mods, wts, grid)
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("time", "surv"))
  expect_equal(nrow(out), length(grid))
})

test_that("blend_survival survival starts at 1 (t = 0)", {
  out <- blend_survival(mods, wts, grid)
  expect_equal(out$surv[out$time == 0], 1, tolerance = 1e-4)
})

test_that("blend_survival survival probabilities lie in [0, 1]", {
  out <- blend_survival(mods, wts, grid)
  expect_true(all(out$surv >= 0 - 1e-8 & out$surv <= 1 + 1e-8))
})

test_that("blend_survival is monotone non-increasing (approx)", {
  out <- blend_survival(mods, wts, grid)
  expect_true(all(diff(out$surv) <= 1e-6))
})

test_that("blend_survival errors on unknown model names in weights_tbl", {
  bad_wts <- wts
  bad_wts$model[1] <- "xdistXXX"
  expect_error(blend_survival(mods, bad_wts, grid), "not found")
})

test_that("predict_surv returns correct length vector", {
  m <- mods[[1]]
  s <- predict_surv(m, times = grid)
  expect_length(s, length(grid))
  expect_type(s, "double")
})
