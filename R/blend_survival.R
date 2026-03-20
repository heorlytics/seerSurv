#' Predict Survival from a Single \code{flexsurvreg} Model
#'
#' @description
#' A thin wrapper around \code{\link[flexsurv]{summary.flexsurvreg}} that
#' extracts the estimated survival probability at each requested time point.
#' Intended for internal use and for users who want individual model curves
#' before blending.
#'
#' @param model A fitted \code{\link[flexsurv]{flexsurvreg}} object.
#' @param times Numeric vector of time points at which to evaluate survival
#'   (same units as the data used to fit \code{model}).
#'
#' @return Numeric vector of survival estimates, one per element of \code{times}.
#'
#' @seealso \code{\link{blend_survival}}
#'
#' @examples
#' \donttest{
#' t <- 0:5
#' s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
#' ipd  <- prep_ipd(t, s)
#' mods <- fit_models(ipd[-1, ])
#' predict_surv(mods[["weibull"]], times = seq(0, 5, by = 0.5))
#' }
#'
#' @export
predict_surv <- function(model, times) {

  if (!inherits(model, "flexsurvreg")) {
    rlang::abort("`model` must be a `flexsurvreg` object.")
  }
  if (!is.numeric(times) || any(times < 0)) {
    rlang::abort("`times` must be a non-negative numeric vector.")
  }

  summary(model, t = times, type = "survival")[[1L]]$est
}


#' Blend Parametric Survival Curves Using Model-Averaging Weights
#'
#' @description
#' Produces a single, weighted-average survival curve by combining predicted
#' survival probabilities from the top-\eqn{k} parametric models selected by
#' \code{\link{compute_weights}}.  The blend is a convex combination:
#' \deqn{\hat{S}(t) = \sum_{i=1}^{k} w_i \, S_i(t),}
#' where \eqn{S_i(t)} is the predicted survival from the \eqn{i}-th model and
#' \eqn{w_i} are the relative-likelihood weights.
#'
#' @param models A named list of fitted \code{\link[flexsurv]{flexsurvreg}}
#'   objects, typically from \code{\link{fit_models}}.
#' @param weights_tbl A tibble from \code{\link{compute_weights}} with columns
#'   \code{model} (matching names in \code{models}) and \code{weight}.
#' @param times Numeric vector of time points at which to evaluate the blended
#'   survival curve.
#'
#' @return A \code{\link[tibble]{tibble}} with two columns:
#' \describe{
#'   \item{\code{time}}{Time points (same as \code{times}).}
#'   \item{\code{surv}}{Blended survival probability.}
#' }
#'
#' @examples
#' \donttest{
#' t <- 0:5
#' s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
#' ipd  <- prep_ipd(t, s)
#' mods <- fit_models(ipd[-1, ])
#' ic   <- extract_ic(mods, "AIC")
#' wts  <- compute_weights(ic, top_k = 3, criterion = "AIC")
#' grid <- seq(0, 10, by = 0.25)
#' S_blend <- blend_survival(mods, wts, grid)
#' head(S_blend)
#' }
#'
#' @export
blend_survival <- function(models, weights_tbl, times) {

  # ---- input validation ---------------------------------------------------- #
  if (!is.list(models) || is.null(names(models))) {
    rlang::abort("`models` must be a named list of flexsurvreg objects.")
  }
  if (!is.data.frame(weights_tbl) ||
      !all(c("model", "weight") %in% names(weights_tbl))) {
    rlang::abort(
      "`weights_tbl` must be a data frame with columns 'model' and 'weight'."
    )
  }
  unknown <- setdiff(weights_tbl$model, names(models))
  if (length(unknown) > 0) {
    rlang::abort(
      paste0("Models named in `weights_tbl` not found in `models`: ",
             paste(unknown, collapse = ", "))
    )
  }
  if (!is.numeric(times) || any(times < 0)) {
    rlang::abort("`times` must be a non-negative numeric vector.")
  }
  # -------------------------------------------------------------------------- #

  models_use <- models[weights_tbl$model]

  surv_mat <- purrr::map_dfc(
    models_use,
    ~ tibble::tibble(surv = predict_surv(.x, times))
  )

  blended_surv <- as.matrix(surv_mat) %*% weights_tbl$weight

  tibble::tibble(
    time = times,
    surv = as.numeric(blended_surv)
  )
}
