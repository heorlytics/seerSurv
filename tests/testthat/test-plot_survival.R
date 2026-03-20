## ============================================================================
## tests/testthat/test-plot_survival.R
## ============================================================================

skip_if_not_installed("ggplot2")

# Minimal fixtures ----------------------------------------------------------- #
t_grid <- seq(0, 5, by = 1)
S_lr   <- tibble::tibble(time = t_grid, surv = c(1, 0.99, 0.98, 0.97, 0.96, 0.95))
S_dr   <- tibble::tibble(time = t_grid, surv = c(1, 0.55, 0.44, 0.39, 0.36, 0.34))
# ---------------------------------------------------------------------------- #

test_that("plot_survival_curves returns a ggplot object", {
  p <- plot_survival_curves(S_lr, S_dr)
  expect_s3_class(p, "ggplot")
})

test_that("plot_survival_curves accepts NULL net curves", {
  p <- plot_survival_curves(S_lr, S_dr, U_L = NULL, U_D = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_survival_curves accepts net curves", {
  U_lr <- tibble::tibble(time = t_grid, surv = S_lr$surv * 0.95)
  U_dr <- tibble::tibble(time = t_grid, surv = S_dr$surv * 0.95)
  p    <- plot_survival_curves(S_lr, S_dr, U_lr, U_dr, tumour = "Test")
  expect_s3_class(p, "ggplot")
})

test_that("plot_survival_curves uses tumour label in title", {
  p <- plot_survival_curves(S_lr, S_dr, tumour = "Melanoma")
  expect_true(grepl("Melanoma", p$labels$title))
})
