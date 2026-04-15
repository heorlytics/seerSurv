predict_surv <- function(model, times) {
  summary(
    model,
    t = times,
    type = "survival"
  )[[1]]$est
}
blend_survival <- function(models, weights_tbl, times) {
  models_use <- models[weights_tbl$model]
  surv_mat <- purrr::map_dfc(
    models_use,
    ~ tibble(surv = predict_surv(.x, times))
  )
  blended_surv <- as.matrix(surv_mat) %*% weights_tbl$weight
  tibble(
    time = times,
    surv = as.numeric(blended_surv)
  )
}
compute_weights <- function(ic_vals, top_k = 3, criterion = "BIC") {

  ranked    <- sort(ic_vals)
  shortlist <- head(ranked, top_k)

  ic_min <- min(shortlist)
  R      <- exp((ic_min - shortlist) / 2)
  w      <- R / sum(R)

  tibble(
    model   = names(shortlist),
    weight  = as.numeric(w),
    ic      = as.numeric(shortlist),
    criterion = criterion
  )
}
extract_ic <- function(models, criterion = "AIC") {
  if (criterion == "AIC") {
    purrr::map_dbl(models, ~ .x$AIC)
  } else {
    purrr::map_dbl(models, ~ .x$BIC)
  }
}
fit_models <- function(ipd) {
  setNames(
    lapply(dists, function(d)
      flexsurvreg(
        Surv(time, event) ~ 1,
        data = ipd,
        weights = weight,
        dist = d
      )
    ),
    dists
  )
}
make_background_surv <- function(lifetable, mean_age, rate_col,time_grid_years) {
  lt <- lifetable %>% filter(Age >= floor(mean_age))

  raw <- tibble(
    time = c(0, seq_len(nrow(lt))),
    surv = c(1, cumprod(1 - lt[[rate_col]]))
  )

  tibble(
    time = time_grid_years,
    surv = approx(raw$time, raw$surv, time_grid_years, rule = 2)$y
  )
}
rmst <- function(time, surv) {
  sum(diff(time) * (head(surv, -1) + tail(surv, -1)) / 2)
}


##wrapper


run_tumour_analysis <- function(
    tumour,
    surv_vec_LR,
    surv_vec_DR,
    N_L,
    N_D,
    prop_male_LR,
    mean_age_LR,
    prop_male_DR,
    mean_age_DR,
    lifetable,
    scenario,
    criterion,
    top_k = 3,
    max_age = 100
) {

  horizon_years <- switch(
    scenario,
    lifetime = max_age - min(mean_age_LR, mean_age_DR),
    `20y`     = 20,
    `5y`      = 5
  )
  time_years <- 0:(length(surv_vec_LR) - 1)
  time_grid_years <- seq(0, horizon_years, by = 1)

  prep_ipd <- function(time, Su) {
    w <- abs(Su - dplyr::lead(Su))
    w[length(w)] <- 1 - sum(w[-length(w)])

    tibble(
      time   = time,
      weight = w,
      event  = c(rep(1, length(time) - 1), 0)
    ) %>% filter(weight > 0)
  }

  ipd_L <- prep_ipd(time_years, surv_vec_LR)
  ipd_D <- prep_ipd(time_years, surv_vec_DR)

  lifetable <- lifetable %>%
    mutate(
      btrate_LR = Males * prop_male_LR + Females * (1 - prop_male_LR),
      btrate_DR = Males * prop_male_DR + Females * (1 - prop_male_DR)
    )

  B_LR <- make_background_surv(lifetable, mean_age_LR, "btrate_LR", time_grid_years)
  B_DR <- make_background_surv(lifetable, mean_age_DR, "btrate_DR", time_grid_years)

  models_L <- fit_models(ipd_L[-1, ])
  models_D <- fit_models(ipd_D[-1, ])

  use_ic <- criterion

  ic_L <- extract_ic(models_L, use_ic)
  ic_D <- extract_ic(models_D, use_ic)

  weights_L <- compute_weights(ic_L, top_k, use_ic)
  weights_D <- compute_weights(ic_D, top_k, use_ic)

  S_blend_L <- blend_survival(models_L, weights_L, time_grid_years)
  S_blend_D <- blend_survival(models_D, weights_D, time_grid_years)

  U_L <- S_blend_L %>%
    left_join(B_LR, by = "time") %>%
    mutate(surv = surv.x * surv.y) %>%
    select(time, surv)

  U_D <- S_blend_D %>%
    left_join(B_DR, by = "time") %>%
    mutate(surv = surv.x * surv.y) %>%
    select(time, surv)

  area_surv <- function(x, y) {
    step <- x[2] - x[1]
    step * (sum(y[-c(1, length(y))]) + 0.5 * (y[1] + y[length(y)]))
  }

  RMST_L_years <- area_surv(U_L$time, U_L$surv)
  RMST_D_years <- area_surv(U_D$time, U_D$surv)

  M_years  <- RMST_L_years - RMST_D_years
  M_months <- M_years * 12

  N_months <- max(time_grid_years) * 12

  objective_PLL <- function(p, M, N) {
    ((1 - p^(N + 1)) / (1 - p) - M)^2
  }

  P_LL <- optimize(
    f = objective_PLL,
    interval = c(0, 0.9999),
    M = M_months,
    N = N_months
  )$minimum

  months <- 0:N_months

  Q_L_m <- approx(time_grid_years * 12, S_blend_L$surv, xout = months, rule = 2)$y
  Q_D_m <- approx(time_grid_years * 12, S_blend_D$surv, xout = months, rule = 2)$y

  ## ---- Approach 1 ----
  t_star <- if (min(Q_L_m) > 0.5) NA else which(Q_L_m <= 0.5)[1]
  L_t <- function(t) P_LL^t

  if (is.na(t_star)) {
    cat("Median not reached in criteria 1")
    P_DL_1 <- NA_real_
    P_Death_L_1 <- NA_real_

  } else {

    denom_A1 <- sum(
      sapply(0:(t_star - 1), function(t)
        L_t(t) * Q_D_m[t_star - t])
    ) + L_t(t_star)

    P_DL_1 <- (0.5 - L_t(t_star)) / denom_A1
    P_Death_L_1 <- 1 - P_LL - P_DL_1
  }

  ## ---- Approach 2 ----
  D_t <- function(t, P_DL) {
    if (t == 0) return(0)
    sum(sapply(0:(t - 1), function(i) L_t(i) * Q_D_m[t - i])) * P_DL
  }

  R_t <- approx(
    time_grid_years * 12,
    S_blend_L$surv,
    xout = months,
    rule = 2
  )$y

  objective_PDL <- function(P_DL) {
    ssq <- 0
    for (t in 1:N_months) {
      pred <- L_t(t) + D_t(t, P_DL)
      obs  <- R_t[t + 1]
      if (obs > 0) ssq <- ssq + (1 - pred / obs)^2
    }
    ssq
  }

  P_DL_2 <- optimize(
    f = objective_PDL,
    interval = c(0, 1 - P_LL)
  )$minimum

  P_Death_L_2 <- 1 - P_LL - P_DL_2

  summary_text <- sprintf(
    paste(
      "%s (LR mean age %.0f): lifetime horizon %.0f years.",
      "Mean sojourn time %.1f months (P(LR)=%.3f).",
      "P(DR)=%.3f (Approach 1) and %.3f (Approach 2)."
    ),
    tumour,
    mean_age_LR,
    max(time_grid_years),
    M_months,
    P_LL,
    P_DL_1,
    P_DL_2
  )

  tibble(
    P_LL = P_LL,
    # P_DL_approach1 = P_DL_1,
    P_DL_approach2 = P_DL_2,
    # P_Death_L_approach1 = P_Death_L_1,
    P_Death_L_approach2 = P_Death_L_2,
    RMST_L_years = RMST_L_years,
    RMST_D_years = RMST_D_years,
    M_years = M_years,
    M_months = M_months,
    # summary_text = summary_text
  )
}
