#' Compute the Restricted Mean Survival Time (Trapezoidal Rule)
#'
#' @description
#' Calculates the area under a survival curve up to the last observed time
#' point using the composite trapezoidal rule.  The result equals the
#' restricted mean survival time (RMST) when the input curve has been
#' evaluated on a sufficiently fine grid.
#'
#' @param time Numeric vector of time points (strictly increasing, same units
#'   as the survival data).
#' @param surv Numeric vector of survival probabilities evaluated at each
#'   element of \code{time}.  Must be the same length as \code{time}.
#'
#' @return A single numeric value: the area under the curve (RMST) in the
#'   same units as \code{time}.
#'
#' @details
#' The trapezoidal approximation is:
#' \deqn{\text{RMST} \approx \sum_{i=1}^{n-1} (t_{i+1} - t_i)
#'   \cdot \frac{S(t_i) + S(t_{i+1})}{2}.}
#' For uniform grids this simplifies to:
#' \deqn{\Delta t \cdot \Bigl[\tfrac{1}{2}S(t_0) + S(t_1) + \cdots +
#'   S(t_{n-2}) + \tfrac{1}{2}S(t_{n-1})\Bigr].}
#'
#' @seealso \code{\link{run_tumour_analysis}}
#'
#' @examples
#' # RMST from a synthetic exponential curve
#' t <- seq(0, 10, by = 0.1)
#' s <- exp(-0.1 * t)
#' rmst(t, s)   # close to 1/0.1 = 10 as horizon -> infinity
#'
#' @export
rmst <- function(time, surv) {

  # ---- input validation ---------------------------------------------------- #
  if (!is.numeric(time) || !is.numeric(surv)) {
    rlang::abort("`time` and `surv` must be numeric vectors.")
  }
  if (length(time) != length(surv)) {
    rlang::abort("`time` and `surv` must have the same length.")
  }
  if (length(time) < 2L) {
    rlang::abort("At least two time points are required.")
  }
  # -------------------------------------------------------------------------- #

  sum(diff(time) * (utils::head(surv, -1L) + utils::tail(surv, -1L)) / 2)
}
