#' Prepare Pseudo-Individual Patient Data from Aggregate Survival Proportions
#'
#' @description
#' Converts a vector of aggregate survival proportions (typically read from a
#' Kaplan–Meier table or a registry publication) into a pseudo-individual
#' patient data (pseudo-IPD) tibble suitable for passing to
#' \code{\link{fit_models}}.  Each row represents one "pseudo-patient"; its
#' weight is the probability mass that the Kaplan–Meier curve drops at that
#' time point.  The last observed subject is treated as censored.
#'
#' @param time Numeric vector of observed time points (in years).  Must be
#'   non-negative, strictly increasing, and the same length as \code{surv}.
#' @param surv Numeric vector of survival probabilities at each time point.
#'   Must start at 1 (or very close to it), be non-increasing, and lie in
#'   \eqn{[0, 1]}.
#'
#' @return A \code{\link[tibble]{tibble}} with three columns:
#' \describe{
#'   \item{\code{time}}{Time point (same units as input).}
#'   \item{\code{weight}}{Probability mass at that time point (non-negative,
#'     sums to 1).}
#'   \item{\code{event}}{Event indicator: \code{1} for an observed event,
#'     \code{0} for the last (censored) observation.}
#' }
#' Rows with zero weight are silently dropped.
#'
#' @details
#' The weight for time \eqn{t_i} is computed as
#' \deqn{w_i = | S(t_i) - S(t_{i+1}) |,}
#' with the final weight adjusted so that all weights sum to 1.  This ensures
#' that the empirical distribution implied by the pseudo-IPD exactly reproduces
#' the input survival curve.
#'
#' @examples
#' # Five-year melanoma LR survival (SEER-17)
#' t <- 0:5
#' s <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
#' ipd <- prep_ipd(t, s)
#' ipd
#'
#' @export
prep_ipd <- function(time, surv) {

  # ---- input validation ---------------------------------------------------- #
  if (!is.numeric(time) || !is.numeric(surv)) {
    rlang::abort("`time` and `surv` must both be numeric vectors.")
  }
  if (length(time) != length(surv)) {
    rlang::abort("`time` and `surv` must have the same length.")
  }
  if (any(diff(time) <= 0)) {
    rlang::abort("`time` must be strictly increasing.")
  }
  if (any(surv < 0 | surv > 1 + .Machine$double.eps)) {
    rlang::abort("`surv` must lie within [0, 1].")
  }
  if (any(diff(surv) > .Machine$double.eps * 100)) {
    rlang::warn("`surv` does not appear to be monotone non-increasing.")
  }
  # -------------------------------------------------------------------------- #

  n <- length(time)
  w <- abs(surv - dplyr::lead(surv))

  # Final weight absorbs any residual so that sum(w) == 1
  w[n] <- max(0, 1 - sum(w[-n], na.rm = TRUE))

  tibble::tibble(
    time   = time,
    weight = w,
    event  = c(rep(1L, n - 1L), 0L)
  ) |>
    dplyr::filter(.data$weight > 0)
}
