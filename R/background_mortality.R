#' Construct Background (All-Cause) Survival from a US-CDC Life-Table
#'
#' @description
#' Uses age- and sex-specific annual death rates from a life-table to build a
#' cumulative background survival curve, then interpolates it onto a user-
#' supplied time grid.  The resulting curve is multiplied against a relative
#' (disease-specific) survival curve inside \code{\link{run_tumour_analysis}}
#' to obtain net (all-cause) survival.
#'
#' @param lifetable A data frame containing the life-table.  Must include an
#'   \code{Age} column (integer ages 0, 1, 2, …) and at least one numeric
#'   column whose name is supplied via \code{rate_col}.  The bundled
#'   \code{\link{lifetable_seer}} dataset is the canonical input.
#' @param mean_age Numeric scalar.  The mean (or median) age of the study
#'   cohort.  The life-table is sub-set to rows where
#'   \code{Age >= floor(mean_age)}.
#' @param rate_col Character scalar.  Name of the column in \code{lifetable}
#'   that contains the annual death rates to use (e.g. \code{"btrate_LR"} or
#'   \code{"btrate_DR"} after sex-mixing has been applied).
#' @param time_grid_years Numeric vector of time points (in years) onto which
#'   the cumulative survival is interpolated.
#'
#' @return A \code{\link[tibble]{tibble}} with two columns:
#' \describe{
#'   \item{\code{time}}{Time in years (same as \code{time_grid_years}).}
#'   \item{\code{surv}}{Background survival probability at each time point.}
#' }
#'
#' @details
#' The cumulative background survival at integer year \eqn{t} (starting from
#' \code{floor(mean_age)}) is:
#' \deqn{B(t) = \prod_{j=0}^{t-1} \bigl(1 - q_j\bigr),}
#' where \eqn{q_j} is the annual death rate at age \eqn{j + \lfloor\text{mean age}\rfloor}.
#' Non-integer time points are linearly interpolated; times beyond the last
#' life-table row use flat extrapolation (\code{rule = 2} in
#' \code{\link[stats]{approx}}).
#'
#' @seealso \code{\link{lifetable_seer}}, \code{\link{run_tumour_analysis}}
#'
#' @examples
#' data(lifetable_seer)
#'
#' # Add a sex-mixed rate column (60% male cohort)
#' lt <- lifetable_seer
#' lt$btrate_mix <- lt$Males * 0.60 + lt$Females * 0.40
#'
#' grid <- seq(0, 39, by = 1)   # 39-year horizon from age 61
#' B <- make_background_surv(lt, mean_age = 61, rate_col = "btrate_mix",
#'                           time_grid_years = grid)
#' head(B)
#'
#' @export
make_background_surv <- function(lifetable, mean_age, rate_col,
                                 time_grid_years) {

  # ---- input validation ---------------------------------------------------- #
  if (!is.data.frame(lifetable)) {
    rlang::abort("`lifetable` must be a data frame.")
  }
  if (!("Age" %in% names(lifetable))) {
    rlang::abort("`lifetable` must contain an 'Age' column.")
  }
  if (!(rate_col %in% names(lifetable))) {
    rlang::abort(
      paste0("`rate_col` '", rate_col, "' not found in `lifetable`.")
    )
  }
  if (!is.numeric(mean_age) || length(mean_age) != 1L || mean_age < 0) {
    rlang::abort("`mean_age` must be a single non-negative number.")
  }
  if (!is.numeric(time_grid_years) || any(time_grid_years < 0)) {
    rlang::abort("`time_grid_years` must be a non-negative numeric vector.")
  }
  # -------------------------------------------------------------------------- #

  lt <- dplyr::filter(lifetable, .data$Age >= floor(mean_age))

  raw <- tibble::tibble(
    time = c(0, seq_len(nrow(lt))),
    surv = c(1, cumprod(1 - lt[[rate_col]]))
  )

  tibble::tibble(
    time = time_grid_years,
    surv = stats::approx(raw$time, raw$surv, time_grid_years, rule = 2L)$y
  )
}
