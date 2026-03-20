## ============================================================================
## tests/testthat/test-background_mortality.R
## ============================================================================

lt <- lifetable_seer
lt$btrate_mix <- lt$Males * 0.58 + lt$Females * 0.42
grid <- seq(0, 39, by = 1)

test_that("make_background_surv returns tibble with time and surv", {
  out <- make_background_surv(lt, 61, "btrate_mix", grid)
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("time", "surv"))
  expect_equal(nrow(out), length(grid))
})

test_that("make_background_surv starts at 1 at time 0", {
  out <- make_background_surv(lt, 61, "btrate_mix", grid)
  expect_equal(out$surv[out$time == 0], 1, tolerance = 1e-6)
})

test_that("make_background_surv survival is non-increasing", {
  out <- make_background_surv(lt, 61, "btrate_mix", grid)
  expect_true(all(diff(out$surv) <= 1e-8))
})

test_that("make_background_surv survival is within (0, 1]", {
  out <- make_background_surv(lt, 61, "btrate_mix", grid)
  expect_true(all(out$surv > 0 & out$surv <= 1))
})

test_that("make_background_surv errors when rate_col missing", {
  expect_error(
    make_background_surv(lifetable_seer, 61, "nonexistent_col", grid),
    "not found"
  )
})

test_that("make_background_surv errors when Age column missing", {
  bad_lt <- lifetable_seer
  names(bad_lt)[names(bad_lt) == "Age"] <- "age_lower"
  expect_error(
    make_background_surv(bad_lt, 61, "Males", grid),
    "Age"
  )
})
