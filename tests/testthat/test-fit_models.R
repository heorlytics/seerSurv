## ============================================================================
## tests/testthat/test-fit_models.R
## ============================================================================

skip_if_not_installed("flexsurv")

t <- 0:5
s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd <- prep_ipd(t, s)

test_that("fit_models returns a named list", {
  mods <- fit_models(ipd[-1, ])
  expect_type(mods, "list")
  expect_named(mods)
})

test_that("fit_models returns flexsurvreg objects", {
  mods <- fit_models(ipd[-1, ])
  expect_true(all(vapply(mods, inherits, logical(1L), "flexsurvreg")))
})

test_that("fit_models returns at most length(dists) models", {
  dists <- c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")
  mods  <- fit_models(ipd[-1, ], dists = dists)
  expect_lte(length(mods), length(dists))
})

test_that("fit_models subset of dists works", {
  mods <- fit_models(ipd[-1, ], dists = c("exp", "weibull"))
  expect_lte(length(mods), 2L)
  expect_true(all(names(mods) %in% c("exp", "weibull")))
})

test_that("fit_models errors on non-data-frame input", {
  expect_error(fit_models(list(time = 1:3, event = c(1, 1, 0), weight = c(0.3, 0.4, 0.3))),
               "data frame")
})

test_that("fit_models errors on missing required columns", {
  bad_ipd <- ipd[-1, ]
  bad_ipd$weight <- NULL
  expect_error(fit_models(bad_ipd), "weight")
})

test_that("fit_models errors when time = 0 rows present", {
  expect_error(fit_models(ipd), "strictly positive")
})

test_that("fit_models errors on empty dists vector", {
  expect_error(fit_models(ipd[-1, ], dists = character(0)), "non-empty")
})
