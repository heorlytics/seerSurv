#' Extract Information Criteria from a List of Parametric Survival Models
#'
#' @description
#' Extracts either the Akaike Information Criterion (AIC) or the Bayesian
#' Information Criterion (BIC) from each \code{flexsurvreg} object in a named
#' list, returning a named numeric vector for use with
#' \code{\link{compute_weights}}.
#'
#' @param models A named list of fitted \code{\link[flexsurv]{flexsurvreg}}
#'   objects, typically the output of \code{\link{fit_models}}.
#' @param criterion Character scalar: \code{"AIC"} (default) or \code{"BIC"}.
#'
#' @return A named numeric vector of information criterion values, one per
#'   model.
#'
#' @examples
#' \donttest{
#' t <- 0:5
#' s <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
#' ipd   <- prep_ipd(t, s)
#' mods  <- fit_models(ipd[-1, ])
#' extract_ic(mods, "AIC")
#' extract_ic(mods, "BIC")
#' }
#'
#' @export
extract_ic <- function(models, criterion = c("AIC", "BIC")) {

  criterion <- match.arg(criterion)

  if (!is.list(models) || length(models) == 0) {
    rlang::abort("`models` must be a non-empty named list of flexsurvreg objects.")
  }

  if (criterion == "AIC") {
    purrr::map_dbl(models, ~ .x$AIC)
  } else {
    purrr::map_dbl(models, ~ .x$BIC)
  }
}


#' Compute Relative-Likelihood Weights for Model Averaging
#'
#' @description
#' Selects the top-\eqn{k} models by information criterion value and assigns
#' relative-likelihood weights using the Akaike / Schwarz weight formula:
#' \deqn{w_i = \frac{\exp\!\bigl((\mathrm{IC}_{\min} - \mathrm{IC}_i)/2\bigr)}
#'              {\sum_{j=1}^{k} \exp\!\bigl((\mathrm{IC}_{\min} - \mathrm{IC}_j)/2\bigr)}.}
#'
#' @param ic_vals Named numeric vector of information criterion values, e.g.
#'   from \code{\link{extract_ic}}.
#' @param top_k Integer.  Number of top-ranked models to retain (default 3).
#' @param criterion Character scalar passed through for labelling only.
#'   Typically \code{"AIC"} or \code{"BIC"}.
#'
#' @return A \code{\link[tibble]{tibble}} with four columns:
#' \describe{
#'   \item{\code{model}}{Distribution name (from names of \code{ic_vals}).}
#'   \item{\code{weight}}{Relative-likelihood model weight (sums to 1).}
#'   \item{\code{ic}}{Raw information criterion value.}
#'   \item{\code{criterion}}{Label passed via the \code{criterion} argument.}
#' }
#'
#' @examples
#' ic <- c(exp = 240.1, weibull = 238.7, gompertz = 241.5,
#'         llogis = 239.3, lnorm = 243.0, gamma = 238.9)
#' compute_weights(ic, top_k = 3, criterion = "AIC")
#'
#' @export
compute_weights <- function(ic_vals, top_k = 3L, criterion = "AIC") {

  # ---- input validation ---------------------------------------------------- #
  if (!is.numeric(ic_vals) || is.null(names(ic_vals))) {
    rlang::abort("`ic_vals` must be a named numeric vector.")
  }
  top_k <- as.integer(top_k)
  if (top_k < 1L || top_k > length(ic_vals)) {
    rlang::abort(
      paste0("`top_k` must be between 1 and ", length(ic_vals),
             " (the number of models).")
    )
  }
  # -------------------------------------------------------------------------- #

  ranked    <- sort(ic_vals)
  shortlist <- utils::head(ranked, top_k)

  ic_min <- min(shortlist)
  R      <- exp((ic_min - shortlist) / 2)
  w      <- R / sum(R)

  tibble::tibble(
    model     = names(shortlist),
    weight    = as.numeric(w),
    ic        = as.numeric(shortlist),
    criterion = criterion
  )
}
