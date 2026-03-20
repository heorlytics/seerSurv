## ============================================================================
## tests/testthat/test-rmst.R
## ============================================================================

test_that("rmst matches analytical result for exponential curve", {
  # For S(t) = exp(-lambda * t), RMST over [0, T] = (1 - exp(-lambda*T)) / lambda
  lambda <- 0.1
  T_max  <- 50
  t      <- seq(0, T_max, by = 0.01)
  s      <- exp(-lambda * t)

  analytical <- (1 - exp(-lambda * T_max)) / lambda
  expect_equal(rmst(t, s), analytical, tolerance = 1e-3)
})

test_that("rmst returns a single numeric scalar", {
  t <- seq(0, 5, by = 1)
  s <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
  out <- rmst(t, s)
  expect_length(out, 1L)
  expect_type(out, "double")
})

test_that("rmst is positive for a declining survival curve", {
  t <- seq(0, 10, by = 0.5)
  s <- exp(-0.2 * t)
  expect_gt(rmst(t, s), 0)
})

test_that("rmst errors on mismatched lengths", {
  expect_error(rmst(1:5, c(1, 0.9, 0.8)), "same length")
})

test_that("rmst errors on single-element input", {
  expect_error(rmst(1, 1), "At least two")
})

test_that("rmst errors on non-numeric time", {
  expect_error(rmst(c("a", "b", "c"), 1:3), "numeric")
})
