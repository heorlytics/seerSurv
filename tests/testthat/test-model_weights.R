## ============================================================================
## tests/testthat/test-model_weights.R
## ============================================================================

ic_vec <- c(exp = 245.1, weibull = 238.7, gompertz = 241.5,
            llogis = 239.3, lnorm  = 243.0, gamma   = 238.9)

test_that("compute_weights returns a tibble with correct columns", {
  out <- compute_weights(ic_vec, top_k = 3, criterion = "AIC")
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("model", "weight", "ic", "criterion"))
})

test_that("compute_weights selects top_k rows", {
  out <- compute_weights(ic_vec, top_k = 3, criterion = "AIC")
  expect_equal(nrow(out), 3L)
})

test_that("compute_weights weights sum to 1", {
  out <- compute_weights(ic_vec, top_k = 3, criterion = "AIC")
  expect_equal(sum(out$weight), 1, tolerance = 1e-10)
})

test_that("compute_weights weights are non-negative", {
  out <- compute_weights(ic_vec, top_k = 3, criterion = "AIC")
  expect_true(all(out$weight >= 0))
})

test_that("compute_weights selects lowest IC models", {
  out <- compute_weights(ic_vec, top_k = 3, criterion = "AIC")
  # weibull (238.7), gamma (238.9), llogis (239.3) should be top 3
  expect_setequal(out$model, c("weibull", "gamma", "llogis"))
})

test_that("compute_weights errors on unnamed vector", {
  expect_error(compute_weights(c(240, 238, 242), top_k = 3), "named")
})

test_that("compute_weights errors on top_k > n models", {
  expect_error(compute_weights(ic_vec, top_k = 10), "between 1")
})

test_that("extract_ic returns named numeric vector for AIC", {
  skip_if_not_installed("flexsurv")
  t <- 0:5
  s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
  ipd  <- prep_ipd(t, s)
  mods <- fit_models(ipd[-1, ])

  aic <- extract_ic(mods, "AIC")
  expect_type(aic, "double")
  expect_named(aic)
})

test_that("extract_ic returns named numeric vector for BIC", {
  skip_if_not_installed("flexsurv")
  t <- 0:5
  s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
  ipd  <- prep_ipd(t, s)
  mods <- fit_models(ipd[-1, ])

  bic <- extract_ic(mods, "BIC")
  expect_type(bic, "double")
  expect_true(all(is.finite(bic)))
})
