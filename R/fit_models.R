#' Fit a Panel of Parametric Survival Models to Pseudo-IPD
#'
#' @description
#' Fits each distribution named in \code{dists} to a pseudo-individual patient
#' data (pseudo-IPD) data frame using \code{\link[flexsurv]{flexsurvreg}} with
#' observation weights.  The resulting named list is designed to be passed
#' directly to \code{\link{extract_ic}} and \code{\link{blend_survival}}.
#'
#' @param ipd A tibble / data frame with at least three columns:
#'   \describe{
#'     \item{\code{time}}{Event or censoring time (numeric, > 0).}
#'     \item{\code{event}}{Event indicator (1 = event, 0 = censored).}
#'     \item{\code{weight}}{Observation weight (non-negative numeric).}
#'   }
#'   Typically the output of \code{\link{prep_ipd}} with the first row
#'   (time = 0) removed.
#' @param dists Character vector of distribution names accepted by
#'   \code{\link[flexsurv]{flexsurvreg}}.  Defaults to the six-distribution
#'   set used throughout \pkg{seerSurv}:
#'   \code{c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")}.
#'
#' @return A named list of \code{flexsurvreg} objects, one per distribution.
#'   Names match the corresponding element of \code{dists}.
#'
#' @details
#' The \code{ipd} data frame should have the first row (time = 0) removed
#' before calling this function because \code{flexsurvreg} requires strictly
#' positive event times.  A common idiom is:
#' \preformatted{
#'   ipd  <- prep_ipd(time_vec, surv_vec)
#'   mods <- fit_models(ipd[-1, ])
#' }
#'
#' @seealso \code{\link{prep_ipd}}, \code{\link{extract_ic}},
#'   \code{\link{blend_survival}}
#'
#' @examples
#' \donttest{
#' t <- 0:5
#' s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)   # melanoma DR (SEER-17)
#' ipd   <- prep_ipd(t, s)
#' mods  <- fit_models(ipd[-1, ])
#' names(mods)
#' }
#'
#' @export
fit_models <- function(
    ipd,
    dists = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")
) {

  # ---- input validation ---------------------------------------------------- #
  if (!is.data.frame(ipd)) {
    rlang::abort("`ipd` must be a data frame or tibble.")
  }
  required_cols <- c("time", "event", "weight")
  missing_cols  <- setdiff(required_cols, names(ipd))
  if (length(missing_cols) > 0) {
    rlang::abort(
      paste0("`ipd` is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    )
  }
  if (!is.character(dists) || length(dists) == 0) {
    rlang::abort("`dists` must be a non-empty character vector.")
  }
  if (any(ipd$time <= 0)) {
    rlang::abort(
      paste0("All `time` values in `ipd` must be strictly positive.  ",
             "Did you forget to remove the time = 0 row? ",
             "Try: fit_models(prep_ipd(t, s)[-1, ])")
    )
  }
  # -------------------------------------------------------------------------- #

  mods <- setNames(
    lapply(dists, function(d) {
      tryCatch(
        flexsurv::flexsurvreg(
          survival::Surv(time, event) ~ 1,
          data    = ipd,
          weights = weight,
          dist    = d
        ),
        error = function(e) {
          rlang::warn(
            paste0("Distribution '", d, "' failed to converge: ", e$message,
                   ".  It will be excluded from IC ranking.")
          )
          NULL
        }
      )
    }),
    dists
  )
  Filter(Negate(is.null), mods)
}
