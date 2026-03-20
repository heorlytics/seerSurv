## ============================================================================
## tests/testthat/test-prep_ipd.R
## ============================================================================

test_that("prep_ipd returns a tibble with correct columns", {
  t <- 0:5
  s <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
  out <- prep_ipd(t, s)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("time", "weight", "event"))
})

test_that("prep_ipd weights sum to 1", {
  t <- 0:5
  s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
  out <- prep_ipd(t, s)

  expect_equal(sum(out$weight), 1, tolerance = 1e-8)
})

test_that("prep_ipd last row has event = 0", {
  t <- 0:5
  s <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
  out <- prep_ipd(t, s)

  expect_equal(tail(out$event, 1), 0L)
})

test_that("prep_ipd drops zero-weight rows", {
  # Flat survival for one interval -> zero weight at that step
  t <- 0:4
  s <- c(1, 0.8, 0.8, 0.6, 0.4)   # no drop at t=2
  out <- prep_ipd(t, s)

  expect_true(all(out$weight > 0))
})

test_that("prep_ipd errors on mismatched lengths", {
  expect_error(prep_ipd(0:5, c(1, 0.9, 0.8)), "same length")
})

test_that("prep_ipd errors on non-numeric input", {
  expect_error(prep_ipd(0:3, c("a", "b", "c", "d")), "numeric")
})

test_that("prep_ipd errors on non-monotone time", {
  expect_error(prep_ipd(c(0, 2, 1, 3), c(1, 0.9, 0.8, 0.7)),
               "strictly increasing")
})
