#' Estimate Cancer-Recurrence Transition Probabilities for One Tumour Type
#'
#' @description
#' The primary user-facing function of \pkg{seerSurv}.  Given five-year aggregate
#' survival proportions for local/regional-recurrence (LR) and distant-
#' recurrence (DR) populations together with demographic inputs, this function:
#' \enumerate{
#'   \item Converts aggregate survival proportions to pseudo-IPD.
#'   \item Fits a panel of parametric models and selects the top-\eqn{k} by
#'         information criterion (AIC or BIC).
#'   \item Blends the selected models using relative-likelihood weights.
#'   \item Adjusts for age- and sex-specific background mortality from
#'         US-CDC life-tables.
#'   \item Computes the mean LR sojourn time (RMST differential) and derives
#'         the monthly LR persistence probability \eqn{P(\text{LR})}.
#'   \item Estimates the monthly LR-to-DR transition probability
#'         \eqn{P(\text{DR})} via a convolution optimisation.
#'   \item Returns all three monthly transition probabilities:
#'         \eqn{P(\text{LR})}, \eqn{P(\text{DR})}, and
#'         \eqn{P(\text{Death} | \text{LR})} = \eqn{1 - P(\text{LR}) - P(\text{DR})}.
#' }
#'
#' @param tumour Character scalar.  Tumour label used for identification in
#'   multi-tumour analyses (e.g. \code{"Melanoma"}).
#' @param surv_vec_LR Numeric vector of survival probabilities for the LR
#'   population at years 0, 1, 2, 3, 4, 5 (must start at 1).
#' @param surv_vec_DR Numeric vector of survival probabilities for the DR
#'   population at years 0, 1, 2, 3, 4, 5 (must start at 1).
#' @param N_L Integer.  Number of LR patients in the SEER cohort (informational;
#'   currently stored for provenance and future sample-size adjustments).
#' @param N_D Integer.  Number of DR patients in the SEER cohort (informational).
#' @param prop_male_LR Numeric in \eqn{[0, 1]}.  Proportion male in the LR
#'   population, used to construct a sex-mixed background mortality rate.
#' @param mean_age_LR Numeric.  Mean (or median) age of the LR population.
#' @param prop_male_DR Numeric in \eqn{[0, 1]}.  Proportion male in the DR
#'   population.
#' @param mean_age_DR Numeric.  Mean (or median) age of the DR population.
#' @param lifetable A data frame with columns \code{Age}, \code{Males},
#'   \code{Females} containing annual death rates.  The bundled
#'   \code{\link{lifetable_seer}} dataset is the canonical choice.
#' @param scenario Character scalar controlling the extrapolation horizon.
#'   One of:
#'   \describe{
#'     \item{\code{"lifetime"}}{Lifetime horizon: \eqn{100 - \min(\text{mean age LR}, \text{mean age DR})} years.}
#'     \item{\code{"20y"}}{Fixed 20-year horizon.}
#'     \item{\code{"5y"}}{Fixed 5-year horizon.}
#'   }
#' @param criterion Character scalar: \code{"AIC"} (default) or \code{"BIC"}.
#'   Governs both model selection and weight computation.
#' @param dists Character vector of distribution names passed to
#'   \code{\link{fit_models}}.  Defaults to the six-distribution base set.
#' @param top_k Integer.  Number of top-ranked models to include in the blend
#'   (default 3).
#' @param max_age Integer.  Maximum age assumed for the lifetime horizon
#'   calculation (default 100).
#'
#' @return A one-row \code{\link[tibble]{tibble}} containing:
#' \describe{
#'   \item{\code{P_LL}}{Monthly probability of remaining in LR.}
#'   \item{\code{P_DL_approach2}}{Monthly probability of transitioning from LR
#'     to DR (convolution optimisation, Approach 2).}
#'   \item{\code{P_Death_L_approach2}}{Monthly probability of death from LR
#'     (\eqn{1 - P_{\text{LL}} - P_{\text{DL}}}).}
#'   \item{\code{RMST_L_years}}{Net RMST for the LR population (years).}
#'   \item{\code{RMST_D_years}}{Net RMST for the DR population (years).}
#'   \item{\code{M_years}}{Mean LR sojourn time (RMST differential, years).}
#'   \item{\code{M_months}}{Mean LR sojourn time (months).}
#' }
#'
#' @details
#' \subsection{P(LR) derivation}{
#' The mean LR sojourn time \eqn{M} (in months) is the area between the
#' net LR and net DR survival curves.  \eqn{P(\text{LR})} is the root of:
#' \deqn{\frac{1 - p^{N+1}}{1 - p} = M,}
#' where \eqn{N} is the total number of monthly time steps in the horizon.
#' }
#'
#' \subsection{P(DR) derivation — Approach 2}{
#' The modelled LR survival \eqn{R(t)} is approximated as the sum of patients
#' still in LR (\eqn{L(t) = P_{\text{LL}}^t}) and those who have progressed to
#' DR:
#' \deqn{R(t) \approx L(t) + P_{\text{DL}} \sum_{i=0}^{t-1} L(i)\, Q_D(t-i),}
#' where \eqn{Q_D(\cdot)} is the monthly DR survival.  \eqn{P_{\text{DL}}} is
#' estimated by minimising the relative sum-of-squares discrepancy between
#' \eqn{R(t)} and the right-hand side over all \eqn{t = 1, \ldots, N}.
#' }
#'
#' @seealso \code{\link{prep_ipd}}, \code{\link{fit_models}},
#'   \code{\link{compute_weights}}, \code{\link{blend_survival}},
#'   \code{\link{make_background_surv}}, \code{\link{tumour_data_seer}}
#'
#' @examples
#' \donttest{
#' data(lifetable_seer)
#'
#' result <- run_tumour_analysis(
#'   tumour       = "Melanoma",
#'   surv_vec_LR  = c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754),
#'   surv_vec_DR  = c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410),
#'   N_L          = 64775L,
#'   N_D          = 3121L,
#'   prop_male_LR = 0.580,
#'   mean_age_LR  = 61,
#'   prop_male_DR = 0.692,
#'   mean_age_DR  = 62,
#'   lifetable    = lifetable_seer,
#'   scenario     = "lifetime",
#'   criterion    = "AIC"
#' )
#' result
#' }
#'
#' @export
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
    scenario   = c("lifetime", "20y", "5y"),
    criterion  = c("AIC", "BIC"),
    dists      = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma"),
    top_k      = 3L,
    max_age    = 100L
) {

  # ---- argument matching ---------------------------------------------------- #
  scenario  <- match.arg(scenario)
  criterion <- match.arg(criterion)

  # ---- input validation ----------------------------------------------------- #
  .check_surv_vec <- function(x, label) {
    if (!is.numeric(x) || length(x) < 2L) {
      rlang::abort(paste0("`", label, "` must be a numeric vector of length >= 2."))
    }
    if (abs(x[1] - 1) > 1e-4) {
      rlang::warn(paste0("`", label, "` should start at 1 (S(0) = 1)."))
    }
  }
  .check_surv_vec(surv_vec_LR, "surv_vec_LR")
  .check_surv_vec(surv_vec_DR, "surv_vec_DR")

  for (nm in c("prop_male_LR", "prop_male_DR")) {
    val <- get(nm)
    if (!is.numeric(val) || val < 0 || val > 1) {
      rlang::abort(paste0("`", nm, "` must be numeric in [0, 1]."))
    }
  }
  if (!is.data.frame(lifetable) ||
      !all(c("Age", "Males", "Females") %in% names(lifetable))) {
    rlang::abort("`lifetable` must be a data frame with columns 'Age', 'Males', 'Females'.")
  }
  top_k <- as.integer(top_k)
  # --------------------------------------------------------------------------- #

  # ---- 1. Time grids -------------------------------------------------------- #
  horizon_years <- switch(
    scenario,
    lifetime = max_age - min(mean_age_LR, mean_age_DR),
    `20y`    = 20,
    `5y`     = 5
  )
  time_years      <- seq_len(length(surv_vec_LR)) - 1L   # 0, 1, ..., n-1
  time_grid_years <- seq(0, horizon_years, by = 1)
  N_months        <- as.integer(max(time_grid_years) * 12)
  months          <- 0:N_months

  # ---- 2. Pseudo-IPD -------------------------------------------------------- #
  ipd_L <- prep_ipd(time_years, surv_vec_LR)
  ipd_D <- prep_ipd(time_years, surv_vec_DR)

  # ---- 3. Sex-mixed background mortality ------------------------------------ #
  lt <- lifetable |>
    dplyr::mutate(
      btrate_LR = .data$Males * prop_male_LR + .data$Females * (1 - prop_male_LR),
      btrate_DR = .data$Males * prop_male_DR + .data$Females * (1 - prop_male_DR)
    )

  B_LR <- make_background_surv(lt, mean_age_LR, "btrate_LR", time_grid_years)
  B_DR <- make_background_surv(lt, mean_age_DR, "btrate_DR", time_grid_years)

  # ---- 4. Fit and select parametric models ---------------------------------- #
  models_L <- fit_models(ipd_L[-1L, ], dists = dists)
  models_D <- fit_models(ipd_D[-1L, ], dists = dists)

  ic_L      <- extract_ic(models_L, criterion)
  ic_D      <- extract_ic(models_D, criterion)
  weights_L <- compute_weights(ic_L, top_k, criterion)
  weights_D <- compute_weights(ic_D, top_k, criterion)

  # ---- 5. Blended relative (disease-specific) survival --------------------- #
  S_blend_L <- blend_survival(models_L, weights_L, time_grid_years)
  S_blend_D <- blend_survival(models_D, weights_D, time_grid_years)

  # ---- 6. Net (all-cause) survival ------------------------------------------ #
  U_L <- S_blend_L |>
    dplyr::left_join(B_LR, by = "time") |>
    dplyr::mutate(surv = .data$surv.x * .data$surv.y) |>
    dplyr::select("time", "surv")

  U_D <- S_blend_D |>
    dplyr::left_join(B_DR, by = "time") |>
    dplyr::mutate(surv = .data$surv.x * .data$surv.y) |>
    dplyr::select("time", "surv")

  # ---- 7. RMST and sojourn time --------------------------------------------- #
  RMST_L_years <- .area_trap(U_L$time, U_L$surv)
  RMST_D_years <- .area_trap(U_D$time, U_D$surv)
  M_years      <- RMST_L_years - RMST_D_years
  M_months     <- M_years * 12

  # ---- 8. Solve for P(LR) --------------------------------------------------- #
  objective_PLL <- function(p, M, N) {
    ((1 - p^(N + 1)) / (1 - p) - M)^2
  }

  P_LL <- stats::optimize(
    f        = objective_PLL,
    interval = c(0, 0.9999),
    M        = M_months,
    N        = N_months
  )$minimum

  # ---- 9. Monthly survival curves ------------------------------------------- #
  Q_L_m <- stats::approx(time_grid_years * 12, S_blend_L$surv,
                         xout = months, rule = 2L)$y
  Q_D_m <- stats::approx(time_grid_years * 12, S_blend_D$surv,
                         xout = months, rule = 2L)$y
  R_t   <- stats::approx(time_grid_years * 12, S_blend_L$surv,
                         xout = months, rule = 2L)$y

  # ---- 10. Estimate background death rate d ---------------------------------- #
  bt  <- dplyr::filter(lt, .data$Age >= min(mean_age_LR, mean_age_DR))
  b_t <- bt$btrate_LR[seq_along(time_grid_years)]
  B_t <- stats::approx(
    x    = time_grid_years * 12,
    y    = b_t,
    xout = months,
    rule = 2L
  )$y

  beta_hat <- sum(-log(B_t[2:(N_months + 1)]) / seq_len(N_months))
  d        <- exp(-beta_hat)

  # ---- 11. Solve for P(DR) — Approach 2 ------------------------------------- #
  L_t <- function(t) P_LL^t

  D_t <- function(t, P_DL) {
    if (t == 0L) return(0)
    P_DL * sum(vapply(
      0:(t - 1L),
      function(i) L_t(i) * Q_D_m[t - i],
      numeric(1L)
    ))
  }

  objective_PDL <- function(P_DL, P_LL, Q_D_m, R_t, N) {
    ssq <- 0
    for (t in seq_len(N)) {
      pred <- L_t(t) + D_t(t, P_DL)
      obs  <- R_t[t + 1L]
      if (obs > 0) {
        ssq <- ssq + (1 - pred / obs)^2
      }
    }
    ssq
  }

  upper_bound <- max(0, 1 - d - P_LL)

  opt_PDL <- stats::optimize(
    f         = objective_PDL,
    interval  = c(0, upper_bound),
    P_LL      = P_LL,
    Q_D_m     = Q_D_m,
    R_t       = R_t,
    N         = N_months
  )

  P_DL_2      <- opt_PDL$minimum
  P_Death_L_2 <- 1 - P_LL - P_DL_2

  # ---- 12. Return results ---------------------------------------------------- #
  tibble::tibble(
    P_LL                 = P_LL,
    P_DL_approach2       = P_DL_2,
    P_Death_L_approach2  = P_Death_L_2,
    RMST_L_years         = RMST_L_years,
    RMST_D_years         = RMST_D_years,
    M_years              = M_years,
    M_months             = M_months
  )
}

# ---- Internal helper: trapezoidal area on a uniform grid -------------------- #
# (Exposed via rmst() for users; this internal version avoids package-level
#  function call overhead inside the optimisation loop.)
.area_trap <- function(x, y) {
  step <- x[2L] - x[1L]
  n    <- length(y)
  step * (sum(y[2L:(n - 1L)]) + 0.5 * (y[1L] + y[n]))
}
