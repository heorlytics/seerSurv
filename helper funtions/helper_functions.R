## ============================================================
##  State-Transition Model — Full Pipeline + Poster Plots
##  horizon accepts "lifetime" OR any numeric value 1–100 yrs
##  POSTER EDITION — 300 dpi, base_size = 22
## ============================================================

library(flexsurv)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)

## ── Distribution shortlist ───────────────────────────────────
dists <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma")

## ── Shared poster colours ────────────────────────────────────
.COL_LR   <- "#004C7D"   # dark navy  — LR
.COL_DR   <- "#5aab2b"   # green      — DR
.COL_PRED <- "#E07B00"   # orange     — Approach 2 predicted
.FILL_M   <- "#D6E8F5"   # light blue — sojourn ribbon
.COL_GRID <- "#D0D9E8"   # axis grid lines

## ── Shared ggplot theme (POSTER EDITION — 300 dpi) ───────────
## base_size bumped to 22 so text reads at arm's length;
## all weights, margins, and key widths scaled accordingly.
.poster_theme <- function(base_size = 22) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(family = "Arial", colour = "#203864"),
      plot.title       = element_text(size   = base_size + 4, face = "bold",
                                      colour = "#004C7D",
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(size   = base_size - 4, colour = "#5a6a8a",
                                      margin = margin(b = 12)),
      axis.title       = element_text(size   = base_size - 2, face = "bold"),
      axis.text        = element_text(size   = base_size - 4),
      legend.title     = element_blank(),
      legend.text      = element_text(size   = base_size - 4),
      legend.position  = "bottom",
      legend.key.width = unit(2.0, "cm"),       # wider swatch — readable at distance
      legend.spacing.x = unit(0.7, "cm"),
      panel.grid.major = element_line(colour = .COL_GRID, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(18, 28, 18, 24),
      strip.text       = element_text(size   = base_size - 1, face = "bold",
                                      colour = "#004C7D")
    )
}

## ── Horizon label helper ─────────────────────────────────────
.horizon_label <- function(horizon_years, horizon_input) {
  if (identical(horizon_input, "lifetime"))
    sprintf("Lifetime (%.0f yr)", horizon_years)
  else
    sprintf("%.0f-year horizon", horizon_years)
}


## ── Helper functions ─────────────────────────────────────────

predict_surv <- function(model, times) {
  summary(model, t = times, type = "survival")[[1]]$est
}

blend_survival <- function(models, weights_tbl, times) {
  models_use <- models[weights_tbl$model]
  surv_mat   <- purrr::map_dfc(models_use,
                               ~ tibble(surv = predict_surv(.x, times)))
  blended    <- as.matrix(surv_mat) %*% weights_tbl$weight
  tibble(time = times, surv = as.numeric(blended))
}

compute_weights <- function(ic_vals, top_k = 3, criterion = "AIC") {
  ranked    <- sort(ic_vals)
  shortlist <- head(ranked, top_k)
  ic_min    <- min(shortlist)
  R         <- exp((ic_min - shortlist) / 2)
  w         <- R / sum(R)
  tibble(model     = names(shortlist),
         weight    = as.numeric(w),
         ic        = as.numeric(shortlist),
         criterion = criterion)
}

extract_ic <- function(models, criterion = "AIC") {
  if (criterion == "AIC") purrr::map_dbl(models, ~ .x$AIC)
  else                    purrr::map_dbl(models, ~ .x$BIC)
}

fit_models <- function(ipd) {
  setNames(
    lapply(dists, function(d)
      flexsurvreg(Surv(time, event) ~ 1, data = ipd,
                  weights = weight, dist = d)),
    dists
  )
}

make_background_surv <- function(lifetable, mean_age, rate_col, time_grid_years) {
  lt  <- lifetable %>% filter(Age >= floor(mean_age))
  raw <- tibble(time = c(0, seq_len(nrow(lt))),
                surv = c(1, cumprod(1 - lt[[rate_col]])))
  tibble(time = time_grid_years,
         surv = approx(raw$time, raw$surv, time_grid_years, rule = 2)$y)
}


## ── Internal plot builders ───────────────────────────────────

## ---- Plot 1: Blended Relative Survival Q_L(t) vs Q_D(t) ----
## What it shows  : parametric-extrapolated relative survival for LR and DR,
##                  shaded to distinguish observed from extrapolated region.
## Why it matters : establishes the raw prognosis gap driving all downstream
##                  transition probability estimates.
.plot_relative_survival <- function(S_blend_L, S_blend_D,
                                    horizon_years, horizon_label,
                                    obs_years,
                                    tumour, criterion) {
  obs_cap <- min(obs_years, horizon_years)

  dat <- bind_rows(
    S_blend_L %>% mutate(Group = "Locoregional (LR)  Q_L(t)"),
    S_blend_D %>% mutate(Group = "Distant (DR)  Q_D(t)")
  ) %>% filter(time <= horizon_years)

  ggplot(dat, aes(x = time, y = surv, colour = Group, linetype = Group)) +
    annotate("rect",
             xmin = obs_cap, xmax = horizon_years,
             ymin = 0, ymax = 1,
             fill = "#F0F4FA", alpha = 0.55) +
    annotate("text",
             x     = obs_cap + (horizon_years - obs_cap) * 0.04,
             y     = 0.94,
             label = "Extrapolation zone",
             size  = 6,                              # ~17 pt at 300 dpi
             colour = "#8A9BBE",
             hjust = 0, fontface = "italic") +
    geom_vline(xintercept = obs_cap,
               colour = "#8A9BBE", linetype = "dotted", linewidth = 1.0) +
    geom_line(linewidth = 1.8) +                     # was 1.3
    scale_colour_manual(
      values = c("Locoregional (LR)  Q_L(t)" = .COL_LR,
                 "Distant (DR)  Q_D(t)"      = .COL_DR)
    ) +
    scale_linetype_manual(
      values = c("Locoregional (LR)  Q_L(t)" = "solid",
                 "Distant (DR)  Q_D(t)"      = "dashed")
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty_breaks(n = 7),
                       limits = c(0, horizon_years),
                       expand = c(0, 0)) +
    labs(
      title    = "Blended Relative Survival \u2014 LR vs DR",
      subtitle = sprintf(
        "%s  |  %s  |  Horizon: %s  |  Dotted line = end of observed data",
        tumour, criterion, horizon_label),
      x = "Time (years)",
      y = "Relative survival  [Q(t)]"
    ) +
    .poster_theme()
}


## ---- Plot 2: Unconditional Survival U_L(t), U_D(t) + M ----
## What it shows  : background-adjusted (all-cause) survival curves.
##                  Shaded gap = M, the mean sojourn time in LR state.
## Why it matters : M links registry data directly to P(L|L) via the identity
##                  [1 - P(L|L)^(N+1)] / [1 - P(L|L)] = M.
.plot_unconditional_survival <- function(U_L, U_D,
                                         M_months, P_LL,
                                         horizon_years, horizon_label,
                                         tumour, criterion) {
  wide <- U_L %>%
    rename(surv_L = surv) %>%
    left_join(U_D %>% rename(surv_D = surv), by = "time") %>%
    filter(time <= horizon_years)

  ann_x <- horizon_years * 0.60
  ann_y <- mean(c(approx(wide$time, wide$surv_L, ann_x)$y,
                  approx(wide$time, wide$surv_D, ann_x)$y))

  dat <- bind_rows(
    U_L %>% mutate(Group = "U_L(t)  LR population"),
    U_D %>% mutate(Group = "U_D(t)  DR population")
  ) %>% filter(time <= horizon_years)

  ggplot(dat, aes(x = time, y = surv)) +
    geom_ribbon(data = wide,
                aes(x = time, ymin = surv_D, ymax = surv_L),
                inherit.aes = FALSE,
                fill = .FILL_M, alpha = 0.65) +
    geom_line(aes(colour = Group, linetype = Group), linewidth = 1.8) +  # was 1.3
    annotate("segment",
             x = ann_x, xend = ann_x,
             y    = approx(wide$time, wide$surv_D, ann_x)$y,
             yend = approx(wide$time, wide$surv_L, ann_x)$y,
             colour = "#004C7D", linewidth = 1.2,               # was 0.9
             arrow = arrow(ends = "both", type = "open",
                           length = unit(0.28, "cm"))) +        # was 0.18
    annotate("label",
             x = ann_x + horizon_years * 0.02, y = ann_y,
             label = sprintf("M = %.1f months\nP(L|L) = %.4f", M_months, P_LL),
             size = 5.5, fontface = "bold", colour = "#004C7D", # was 3.6
             fill = "white", label.size = 0.4, hjust = 0) +    # was 0.3
    scale_colour_manual(
      values = c("U_L(t)  LR population" = .COL_LR,
                 "U_D(t)  DR population" = .COL_DR)
    ) +
    scale_linetype_manual(
      values = c("U_L(t)  LR population" = "solid",
                 "U_D(t)  DR population" = "dashed")
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty_breaks(n = 7),
                       limits = c(0, horizon_years),
                       expand = c(0, 0)) +
    labs(
      title    = "Unconditional Survival & Mean LR Sojourn Time (M)",
      subtitle = sprintf(
        "%s  |  %s  |  Horizon: %s  |  Shaded area = M \u2192 solves for P(L|L)",
        tumour, criterion, horizon_label),
      x = "Time (years)",
      y = "Overall survival  [U(t) = B(t) \u00d7 Q(t)]"
    ) +
    .poster_theme()
}


## ---- Plot 3: Approach 2 fit — D(t)+L(t) vs R(t) ----
## What it shows  : transition-based LR survival [D(t)+L(t)] overlaid on the
##                  blended parametric extrapolation [R(t)], plus a residual
##                  panel showing the relative fit error at each month.
## Why it matters : tight overlap + flat residuals validate that the optimised
##                  P(D|L) reproduces the observed LR survival experience.
.plot_approach2_fit <- function(P_LL, P_DL_2,
                                Q_D_m, R_t,
                                N_months, horizon_years, horizon_label,
                                tumour, criterion) {
  L_t <- function(t) P_LL^t

  ## D_t — identical indexing to the wrapper
  D_t <- function(t, P_DL) {
    if (t == 0L) return(0)
    sum(sapply(0L:(t - 1L), function(i) L_t(i) * Q_D_m[t - i])) * P_DL
  }

  plot_months <- min(N_months, 120L)   # cap at 10 yrs for poster
  mo_seq      <- 0L:plot_months

  pred_surv <- vapply(mo_seq,
                      function(t) L_t(t) + D_t(t, P_DL_2),
                      numeric(1L))
  extrap    <- R_t[mo_seq + 1L]
  resid     <- ifelse(extrap > 0, (pred_surv - extrap) / extrap, NA_real_)

  base_dat <- tibble(
    years        = mo_seq / 12,
    Extrapolated = extrap,
    Predicted    = pred_surv,
    residual     = resid
  )

  long_dat <- base_dat %>%
    select(years, Extrapolated, Predicted) %>%
    pivot_longer(-years, names_to = "Series", values_to = "surv") %>%
    mutate(Series = recode(Series,
                           "Extrapolated" = "R(t) \u2014 blended parametric [LR data]",
                           "Predicted"    = "D(t)+L(t) \u2014 transition-based "))

  p_main <- ggplot(long_dat,
                   aes(x = years, y = surv,
                       colour = Series, linetype = Series)) +
    geom_line(linewidth = 1.7) +                               # was 1.25
    scale_colour_manual(
      values = c("R(t) \u2014 blended parametric [LR data]"       = .COL_LR,
                 "D(t)+L(t) \u2014 transition-based " = .COL_PRED)
    ) +
    scale_linetype_manual(
      values = c("R(t) \u2014 blended parametric [LR data]"       = "solid",
                 "D(t)+L(t) \u2014 transition-based " = "longdash")
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty_breaks(n = 6),
                       limits = c(0, plot_months / 12),
                       expand = c(0, 0)) +
    labs(
      title    = "Predicted vs Extrapolated LR Survival",
      subtitle = sprintf(
        "%s  |  %s  |  Horizon: %s  |  P(L|L)=%.4f  P(D|L)=%.4f  P(Death|L)=%.4f",
        tumour, criterion, horizon_label,
        P_LL, P_DL_2, 1 - P_LL - P_DL_2),
      x = NULL, y = "Survival from LR state"
    ) +
    .poster_theme() +
    theme(legend.position = "bottom",
          axis.text.x    = element_blank(),
          axis.ticks.x   = element_blank(),
          plot.margin    = margin(18, 28, 2, 24))              # was margin(10,16,2,14)

  p_resid <- ggplot(base_dat, aes(x = years, y = residual)) +
    geom_hline(yintercept = 0,
               colour = .COL_LR, linewidth = 1.0) +           # was 0.7
    geom_line(colour = .COL_PRED, linewidth = 1.2) +          # was 0.9
    geom_ribbon(aes(ymin = 0, ymax = residual),
                fill = .COL_PRED, alpha = 0.15) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_x_continuous(breaks = pretty_breaks(n = 6),
                       limits = c(0, plot_months / 12),
                       expand = c(0, 0)) +
    labs(x = "Time (years)", y = "Relative\nresidual") +
    .poster_theme() +
    theme(plot.margin   = margin(2, 28, 18, 24),              # was margin(2,16,10,14)
          plot.title    = element_blank(),
          plot.subtitle = element_blank())

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_main / p_resid + patchwork::plot_layout(heights = c(3, 1))
  } else {
    p_main
  }
}


## ---- Plot 4: Model blend weights ----
## What it shows  : the three distributions selected for LR and DR and their
##                  relative AIC/BIC likelihood weights.
## Why it matters : transparency — shows model uncertainty is handled by
##                  averaging over competing parametric forms, not a single
##                  arbitrary choice.
.plot_blend_weights <- function(weights_L, weights_D,
                                horizon_label, tumour, criterion) {
  dist_labels <- c(
    exp      = "Exponential",
    weibull  = "Weibull",
    lnorm    = "Log-normal",
    llogis   = "Log-logistic",
    gompertz = "Gompertz",
    gengamma = "Gen. Gamma"
  )

  bind_rows(
    weights_L %>% mutate(Population = "Locoregional (LR)"),
    weights_D %>% mutate(Population = "Distant (DR)")
  ) %>%
    mutate(
      dist_name = dplyr::recode(model, !!!dist_labels, .default = toupper(model)),
      dist_name = factor(dist_name, levels = rev(sort(unique(dist_name))))
    ) %>%
    ggplot(aes(y = dist_name, x = weight, fill = Population)) +
    geom_col(position = position_dodge(width = 0.72),
             width = 0.62, colour = "white") +
    geom_text(aes(label = sprintf("%.3f", weight)),
              position = position_dodge(width = 0.72),
              hjust = -0.12, size = 5.5, fontface = "bold",   # was 3.4
              colour = "#203864") +
    scale_fill_manual(
      values = c("Locoregional (LR)" = .COL_LR,
                 "Distant (DR)"      = .COL_DR)
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1.15), expand = c(0, 0)) +
    labs(
      title    = "Model Blend Weights \u2014 Top-3 Distributions",
      subtitle = sprintf(
        "%s  |  %s  |  Horizon: %s  |  Weights from %s relative likelihood",
        tumour, criterion, horizon_label, criterion),
      x = "Blend weight", y = NULL
    ) +
    .poster_theme() +
    theme(legend.position = "bottom")
}


## ---- Plot 5: All-in-one — Q_L, U_L, Q_D, U_D in a single panel ----
## What it shows  : all four survival curves together so the reader can see
##                  (a) the raw prognosis gap [Q_L vs Q_D],
##                  (b) the background mortality drag per population
##                      [gap between Q and U for each group],
##                  (c) the sojourn area M [ribbon between U_L and U_D].
## Why it matters : one self-contained figure that tells the full story —
##                  suitable as the centrepiece results plot on the poster.
.plot_survival_overview <- function(S_blend_L, S_blend_D,
                                    U_L, U_D,
                                    M_months, P_LL,
                                    horizon_years, horizon_label,
                                    obs_years,
                                    tumour, criterion,
                                    display_years = NULL,
                                    show_ribbon   = TRUE) {

  x_max   <- if (!is.null(display_years)) min(display_years, horizon_years)
  else horizon_years
  obs_cap <- min(obs_years, x_max)

  ## ── long data (4 curves) ─────────────────────────────────
  dat <- bind_rows(
    S_blend_L %>% filter(time <= x_max) %>%
      mutate(Group = "LR  \u2014  Q_L(t)  [relative]"),
    S_blend_D %>% filter(time <= x_max) %>%
      mutate(Group = "DR  \u2014  Q_D(t)  [relative]"),
    U_L %>% filter(time <= x_max) %>%
      mutate(Group = "LR  \u2014  U_L(t)  [unconditional]"),
    U_D %>% filter(time <= x_max) %>%
      mutate(Group = "DR  \u2014  U_D(t)  [unconditional]")
  )

  dat$Group <- factor(dat$Group, levels = c(
    "LR  \u2014  Q_L(t)  [relative]",
    "LR  \u2014  U_L(t)  [unconditional]",
    "DR  \u2014  Q_D(t)  [relative]",
    "DR  \u2014  U_D(t)  [unconditional]"
  ))

  ## ── colour / linetype / linewidth maps ───────────────────
  col_map <- c(
    "LR  \u2014  Q_L(t)  [relative]"      = "#6FA8D4",
    "LR  \u2014  U_L(t)  [unconditional]"  = "#004C7D",
    "DR  \u2014  Q_D(t)  [relative]"      = "#9BD468",
    "DR  \u2014  U_D(t)  [unconditional]"  = "#3A8A00"
  )
  lty_map <- c(
    "LR  \u2014  Q_L(t)  [relative]"      = "dashed",
    "LR  \u2014  U_L(t)  [unconditional]"  = "solid",
    "DR  \u2014  Q_D(t)  [relative]"      = "dashed",
    "DR  \u2014  U_D(t)  [unconditional]"  = "solid"
  )
  lwd_map <- c(
    "LR  \u2014  Q_L(t)  [relative]"      = 1.1,   # was 0.9
    "LR  \u2014  U_L(t)  [unconditional]"  = 2.0,   # was 1.5
    "DR  \u2014  Q_D(t)  [relative]"      = 1.1,
    "DR  \u2014  U_D(t)  [unconditional]"  = 2.0
  )

  ## ── ribbon data ──────────────────────────────────────────
  ribbon_LR <- S_blend_L %>%
    rename(surv_Q = surv) %>%
    left_join(U_L %>% rename(surv_U = surv), by = "time") %>%
    filter(time <= x_max)

  ribbon_DR <- S_blend_D %>%
    rename(surv_Q = surv) %>%
    left_join(U_D %>% rename(surv_U = surv), by = "time") %>%
    filter(time <= x_max)

  ## ── M annotation position ────────────────────────────────
  wide_U <- U_L %>% rename(surv_L = surv) %>%
    left_join(U_D %>% rename(surv_D = surv), by = "time") %>%
    filter(time <= x_max)

  ann_x <- x_max * 0.62
  ann_y <- mean(c(approx(wide_U$time, wide_U$surv_L, ann_x)$y,
                  approx(wide_U$time, wide_U$surv_D, ann_x)$y))

  ## ── build plot ───────────────────────────────────────────
  p <- ggplot(dat, aes(x = time, y = surv,
                       colour   = Group,
                       linetype = Group)) +
    ## extrapolation zone
    annotate("rect",
             xmin = obs_cap, xmax = x_max, ymin = 0, ymax = 1,
             fill = "#F0F4FA", alpha = 0.50) +
    annotate("text",
             x = obs_cap + (x_max - obs_cap) * 0.04, y = 0.97,
             label = "Extrapolation", size = 5.5,              # was 3.0
             colour = "#8A9BBE", hjust = 0, fontface = "italic") +
    geom_vline(xintercept = obs_cap,
               colour = "#8A9BBE", linetype = "dotted", linewidth = 0.9)  # was 0.6

  if (show_ribbon) {
    p <- p +
      ## LR background gap (between Q_L and U_L)
      geom_ribbon(data = ribbon_LR,
                  aes(x = time, ymin = surv_U, ymax = surv_Q),
                  inherit.aes = FALSE,
                  fill = "#6FA8D4", alpha = 0.16) +
      ## DR background gap (between Q_D and U_D)
      geom_ribbon(data = ribbon_DR,
                  aes(x = time, ymin = surv_U, ymax = surv_Q),
                  inherit.aes = FALSE,
                  fill = "#9BD468", alpha = 0.16) +
      ## sojourn M area (between U_L and U_D)
      geom_ribbon(data = wide_U,
                  aes(x = time, ymin = surv_D, ymax = surv_L),
                  inherit.aes = FALSE,
                  fill = .FILL_M, alpha = 0.55)
  }

  p <- p +
    ## M double-headed arrow
    annotate("segment",
             x = ann_x, xend = ann_x,
             y    = approx(wide_U$time, wide_U$surv_D, ann_x)$y,
             yend = approx(wide_U$time, wide_U$surv_L, ann_x)$y,
             colour = "#004C7D", linewidth = 1.2,              # was 0.9
             arrow = arrow(ends = "both", type = "open",
                           length = unit(0.25, "cm"))) +       # was 0.17
    annotate("label",
             x     = ann_x + x_max * 0.015,
             y     = ann_y,
             label = sprintf("M = %.1f mo\nP(L|L) = %.4f", M_months, P_LL),
             size = 5.5, fontface = "bold", colour = "#004C7D", # was 3.5
             fill = "white", label.size = 0.35, hjust = 0) +   # was 0.25
    ## four curves — linewidth applied per-group via separate geom_line() calls
    ## to avoid the filled-rectangle legend key that scale_discrete_manual("linewidth")
    ## produces in ggplot2 >= 3.4
    geom_line(data = ~ subset(.x, Group == "LR  \u2014  Q_L(t)  [relative]"),
              linewidth = lwd_map[["LR  \u2014  Q_L(t)  [relative]"]]) +
    geom_line(data = ~ subset(.x, Group == "LR  \u2014  U_L(t)  [unconditional]"),
              linewidth = lwd_map[["LR  \u2014  U_L(t)  [unconditional]"]]) +
    geom_line(data = ~ subset(.x, Group == "DR  \u2014  Q_D(t)  [relative]"),
              linewidth = lwd_map[["DR  \u2014  Q_D(t)  [relative]"]]) +
    geom_line(data = ~ subset(.x, Group == "DR  \u2014  U_D(t)  [unconditional]"),
              linewidth = lwd_map[["DR  \u2014  U_D(t)  [unconditional]"]]) +
    scale_colour_manual(values   = col_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty_breaks(n = 7),
                       limits = c(0, x_max), expand = c(0, 0)) +
    labs(
      title    = "Survival Overview: Relative & Unconditional \u2014 LR vs DR",
      subtitle = sprintf(
        paste0("%s  |  %s  |  Horizon: %s\n",
               "Solid = unconditional U(t)=B(t)\u00d7Q(t)  \u00b7  ",
               "Dashed = relative Q(t)  \u00b7  ",
               "Shaded ribbons: navy = LR background gap, ",
               "green = DR background gap, blue = sojourn M"),
        tumour, criterion, horizon_label),
      x = "Time (years)",
      y = "Survival probability"
    ) +
    .poster_theme() +
    theme(
      legend.position  = "bottom",
      legend.key.width = unit(2.2, "cm"),   # was 1.5
      legend.spacing.x = unit(0.6, "cm")   # was 0.5
    )

  p
}


## ── Main wrapper ─────────────────────────────────────────────

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
    horizon   = "lifetime",  # "lifetime"  OR  numeric years 1–100
    criterion = "AIC",
    top_k     = 3L,
    max_age   = 100
) {

  ## ── 1. Validate & resolve horizon ──────────────────────────
  if (identical(horizon, "lifetime")) {
    horizon_years <- max_age - min(mean_age_LR, mean_age_DR)
  } else {
    if (!is.numeric(horizon) || length(horizon) != 1L ||
        horizon < 1 || horizon > 100) {
      stop('`horizon` must be "lifetime" or a single number between 1 and 100.')
    }
    horizon_years <- horizon
  }

  horizon_label   <- .horizon_label(horizon_years, horizon)
  time_years      <- 0:(length(surv_vec_LR) - 1L)
  obs_years       <- max(time_years)
  time_grid_years <- seq(0, horizon_years, by = 1)

  ## ── 2. Pseudo-IPD ──────────────────────────────────────────
  prep_ipd <- function(time, Su) {
    w            <- abs(Su - dplyr::lead(Su))
    w[length(w)] <- 1 - sum(w[-length(w)])
    tibble(time   = time,
           weight = w,
           event  = c(rep(1L, length(time) - 1L), 0L)) %>%
      filter(weight > 0)
  }

  ipd_L <- prep_ipd(time_years, surv_vec_LR)
  ipd_D <- prep_ipd(time_years, surv_vec_DR)

  ## ── 3. Background survival ─────────────────────────────────
  lifetable <- lifetable %>%
    mutate(
      btrate_LR = Males * prop_male_LR + Females * (1 - prop_male_LR),
      btrate_DR = Males * prop_male_DR + Females * (1 - prop_male_DR)
    )

  B_LR <- make_background_surv(lifetable, mean_age_LR, "btrate_LR", time_grid_years)
  B_DR <- make_background_surv(lifetable, mean_age_DR, "btrate_DR", time_grid_years)

  ## ── 4. Parametric fits ─────────────────────────────────────
  models_L <- fit_models(ipd_L[-1L, ])
  models_D <- fit_models(ipd_D[-1L, ])

  ## ── 5. IC weights + blended curves ────────────────────────
  ic_L      <- extract_ic(models_L, criterion)
  ic_D      <- extract_ic(models_D, criterion)

  weights_L <- compute_weights(ic_L, top_k, criterion)
  weights_D <- compute_weights(ic_D, top_k, criterion)

  S_blend_L <- blend_survival(models_L, weights_L, time_grid_years)
  S_blend_D <- blend_survival(models_D, weights_D, time_grid_years)

  ## ── 6. Unconditional survival ──────────────────────────────
  U_L <- S_blend_L %>%
    left_join(B_LR, by = "time") %>%
    mutate(surv = surv.x * surv.y) %>%
    select(time, surv)

  U_D <- S_blend_D %>%
    left_join(B_DR, by = "time") %>%
    mutate(surv = surv.x * surv.y) %>%
    select(time, surv)

  ## ── 7. Mean sojourn time M (trapezoidal rule) ──────────────
  area_trap <- function(x, y) {
    step <- x[2L] - x[1L]
    step * (sum(y[-c(1L, length(y))]) + 0.5 * (y[1L] + y[length(y)]))
  }

  RMST_L_years <- area_trap(U_L$time, U_L$surv)
  RMST_D_years <- area_trap(U_D$time, U_D$surv)
  M_years      <- RMST_L_years - RMST_D_years
  M_months     <- M_years * 12
  N_months     <- as.integer(round(max(time_grid_years) * 12))

  ## ── 8. Solve for P(L|L) ───────────────────────────────────
  objective_PLL <- function(p, M, N) ((1 - p^(N + 1)) / (1 - p) - M)^2

  P_LL <- optimize(objective_PLL,
                   interval = c(1e-6, 0.9999),
                   M = M_months, N = N_months)$minimum

  ## ── 9. Monthly interpolation ──────────────────────────────
  months <- 0L:N_months

  Q_L_m <- approx(time_grid_years * 12, S_blend_L$surv,
                  xout = months, rule = 2)$y
  Q_D_m <- approx(time_grid_years * 12, S_blend_D$surv,
                  xout = months, rule = 2)$y
  R_t   <- approx(time_grid_years * 12, S_blend_L$surv,
                  xout = months, rule = 2)$y

  L_t <- function(t) P_LL^t

  ## ── 10. Approach 1 ────────────────────────────────────────
  t_star <- if (min(Q_L_m) > 0.5) NA_integer_ else which(Q_L_m <= 0.5)[1L]

  if (is.na(t_star)) {
    message(sprintf("[%s | %s | %s] Median not reached — Approach 1 skipped.",
                    tumour, criterion, horizon_label))
    P_DL_1      <- NA_real_
    P_Death_L_1 <- NA_real_
  } else {
    denom_A1    <- sum(sapply(0L:(t_star - 1L),
                              function(t) L_t(t) * Q_D_m[t_star - t])) +
      L_t(t_star)
    P_DL_1      <- (0.5 - L_t(t_star)) / denom_A1
    P_Death_L_1 <- 1 - P_LL - P_DL_1
  }

  ## ── 11. Approach 2 ────────────────────────────────────────
  D_t <- function(t, P_DL) {
    if (t == 0L) return(0)
    sum(sapply(0L:(t - 1L), function(i) L_t(i) * Q_D_m[t - i])) * P_DL
  }

  objective_PDL <- function(P_DL) {
    ssq <- 0
    for (t in 1L:N_months) {
      pred <- L_t(t) + D_t(t, P_DL)
      obs  <- R_t[t + 1L]
      if (obs > 0) ssq <- ssq + (1 - pred / obs)^2
    }
    ssq
  }

  P_DL_2      <- optimize(objective_PDL,
                          interval = c(0, 1 - P_LL))$minimum
  P_Death_L_2 <- 1 - P_LL - P_DL_2

  ## ── 12. Build all five plots ──────────────────────────────
  p1 <- .plot_relative_survival(
    S_blend_L, S_blend_D,
    horizon_years, horizon_label, obs_years,
    tumour, criterion
  )

  p2 <- .plot_unconditional_survival(
    U_L, U_D,
    M_months, P_LL,
    horizon_years, horizon_label,
    tumour, criterion
  )

  p3 <- .plot_approach2_fit(
    P_LL, P_DL_2,
    Q_D_m, R_t,
    N_months, horizon_years, horizon_label,
    tumour, criterion
  )

  p4 <- .plot_blend_weights(
    weights_L, weights_D,
    horizon_label, tumour, criterion
  )

  p5 <- .plot_survival_overview(
    S_blend_L, S_blend_D,
    U_L, U_D,
    M_months, P_LL,
    horizon_years, horizon_label,
    obs_years,
    tumour, criterion
  )

  ## ── 13. Return ───────────────────────────────────────────
  list(
    results = tibble(
      tumour               = tumour,
      horizon_years        = horizon_years,
      horizon_label        = horizon_label,
      criterion            = criterion,
      P_LL                 = P_LL,
      P_DL_approach2       = P_DL_2,
      P_Death_L_approach2  = P_Death_L_2,
      RMST_L_years         = RMST_L_years,
      RMST_D_years         = RMST_D_years,
      M_years              = M_years,
      M_months             = M_months
    ),

    plots = list(
      p1_relative_survival      = p1,
      p2_unconditional_survival = p2,
      p3_approach2_fit          = p3,
      p4_blend_weights          = p4,
      p5_survival_overview      = p5
    ),

    internals = list(
      S_blend_L     = S_blend_L,
      S_blend_D     = S_blend_D,
      U_L           = U_L,
      U_D           = U_D,
      weights_L     = weights_L,
      weights_D     = weights_D,
      Q_D_m         = Q_D_m,
      R_t           = R_t,
      N_months      = N_months,
      horizon_years = horizon_years,
      obs_years     = obs_years
    )
  )
}


## ── Convenience helpers ──────────────────────────────────────

## Formatted console summary
print_results <- function(run) {
  r <- run$results
  cat(sprintf(
    "\n=== %s | %s | %s ===\n  P(L|L)       = %.4f\n  P(D|L)       = %.4f \n  P(Death|L)   = %.4f\n  M            = %.1f months\n  RMST LR      = %.2f yrs\n  RMST DR      = %.2f yrs\n",
    r$tumour, r$criterion, r$horizon_label,
    r$P_LL, r$P_DL_approach2, r$P_Death_L_approach2,
    r$M_months, r$RMST_L_years, r$RMST_D_years
  ))
}

## ── Save all five plots to PNG at poster resolution ──────────
## Default canvas: 16 × 10 in at 300 dpi  ≈  4800 × 3000 px.
## Approach-2 stacked plot gets +3 in of height automatically.
save_run_plots <- function(run,
                           outdir = ".",
                           width  = 16,      # inches  (was 9)
                           height = 10,      # inches  (was 6)
                           dpi    = 300) {

  r      <- run$results
  prefix <- file.path(
    outdir,
    gsub("[ /()]", "_",
         paste(r$tumour, r$criterion,
               gsub(" ", "_", r$horizon_label), sep = "_"))
  )

  for (nm in names(run$plots)) {
    h <- if (grepl("approach2", nm)) height + 3 else height
    ggsave(paste0(prefix, "_", nm, ".png"),
           plot   = run$plots[[nm]],
           width  = width,
           height = h,
           dpi    = dpi,
           bg     = "white")
    message("Saved: ", basename(paste0(prefix, "_", nm, ".png")))
  }
  invisible(NULL)
}

## ── 2×2 combined panel at poster size ───────────────────────
## Pass save_path to also write a PNG in one call.
## Default canvas: 32 × 22 in at 300 dpi — two 16-in panels side-by-side.
combine_plots <- function(run,
                          title      = NULL,
                          save_path  = NULL,   # e.g. "figures/combined.png"
                          width      = 32,     # was not configurable before
                          height     = 22,
                          dpi        = 300) {

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Install patchwork:  install.packages('patchwork')")

  p        <- run$plots
  combined <- (p$p1_relative_survival | p$p2_unconditional_survival) /
    (p$p3_approach2_fit    | p$p4_blend_weights)

  if (!is.null(title))
    combined <- combined + patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(
        size   = 28,           # was 18
        face   = "bold",
        colour = "#004C7D",
        hjust  = 0.5,
        family = "Arial"))
    )

  if (!is.null(save_path)) {
    ggsave(save_path, plot = combined,
           width = width, height = height, dpi = dpi, bg = "white")
    message("Saved combined panel: ", save_path)
  }

  combined
}


## ══════════════════════════════════════════════════════════════
##  EXAMPLE USAGE
## ══════════════════════════════════════════════════════════════
if (FALSE) {

  ## "lifetime" — auto-computed from max_age - min(mean ages)
  run_lt <- run_tumour_analysis(
    tumour = "Melanoma", ...,
    horizon = "lifetime", criterion = "AIC"
  )

  ## Any fixed horizon: 1 to 100 years
  run_20 <- run_tumour_analysis(..., horizon = 20,  criterion = "AIC")
  run_15 <- run_tumour_analysis(..., horizon = 15,  criterion = "AIC")
  run_10 <- run_tumour_analysis(..., horizon = 10,  criterion = "BIC")
  run_5  <- run_tumour_analysis(..., horizon = 5,   criterion = "BIC")
  run_1  <- run_tumour_analysis(..., horizon = 1,   criterion = "AIC")  # stress test

  ## Access results and plots
  print_results(run_lt)
  run_lt$results
  run_lt$plots$p1_relative_survival
  run_lt$plots$p2_unconditional_survival
  run_lt$plots$p3_approach2_fit        # stacked main + residual
  run_lt$plots$p4_blend_weights
  run_lt$plots$p5_survival_overview

  ## 2×2 panel (display only)
  combine_plots(run_lt, title = "Melanoma — AIC | Lifetime")

  ## 2×2 panel + save in one call
  combine_plots(run_lt,
                title     = "Melanoma — AIC | Lifetime",
                save_path = "figures/melanoma_combined.png")

  ## Save all 5 individual PNGs at 300 dpi / 16×10 in
  save_run_plots(run_lt, outdir = "figures/")

  ## Sweep across all horizons and stack results
  runs <- lapply(
    list("lifetime", 30, 20, 15, 10, 5),
    function(h) run_tumour_analysis(..., horizon = h, criterion = "AIC")
  )

  bind_rows(lapply(runs, `[[`, "results")) %>%
    select(horizon_label, P_LL, P_DL_approach2, M_months)
}
