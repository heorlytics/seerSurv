#' seerSurv: Cancer Progression Transition Probabilities from Aggregate Survival Data
#'
#' @description
#' The \pkg{seerSurv} package provides a transparent, reproducible framework for
#' estimating monthly transition probabilities used in state-transition health
#' economic models that include a local/regional recurrence (LR) state.
#'
#' Survival data for LR and distant-recurrence (DR) populations—typically
#' extracted from the SEER-17 registry—are extrapolated using a weighted blend
#' of parametric survival models.  Age- and sex-specific background mortality
#' from US-CDC life-tables is layered on top to produce net (all-cause) survival
#' curves.  The area differential between the LR and DR curves yields the mean
#' LR sojourn time, from which \eqn{P(\text{LR})} is derived.  A convolution
#' optimisation then recovers \eqn{P(\text{DR})}.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item Prepare pseudo-individual patient data from aggregate survival
#'         proportions with \code{\link{prep_ipd}}.
#'   \item Fit a panel of parametric models with \code{\link{fit_models}}.
#'   \item Select and weight the top-\eqn{k} models by AIC or BIC using
#'         \code{\link{extract_ic}} and \code{\link{compute_weights}}.
#'   \item Blend the weighted survival curves with \code{\link{blend_survival}}.
#'   \item Adjust for background mortality with
#'         \code{\link{make_background_surv}}.
#'   \item Run the complete analysis pipeline for one or many tumour types
#'         with \code{\link{run_tumour_analysis}}.
#' }
#'
#' @section Bundled data:
#' \describe{
#'   \item{\code{\link{lifetable_seer}}}{
#'     US-CDC annual death-rate life-table (ages 0–100) with separate
#'     columns for males and females, aligned to the SEER-17 study period.
#'   }
#'   \item{\code{\link{tumour_data_seer}}}{
#'     Reference 5-year aggregate survival proportions (years 0–5) for LR
#'     and DR populations across 11 cancer types derived from SEER-17, together
#'     with sample sizes, sex proportions, and mean ages.
#'   }
#' }
#'
#' @references
#' Mansoori S, Pandey S, Rani R, Singh B, Kurt M (2025).
#' "Deriving Cancer Progression Rates After Local/Regional Recurrence Using
#' Aggregate Survival Data: An R Package for Health Economic Evaluations."
#' \emph{Value in Health} (submitted).
#'
#' Surveillance, Epidemiology, and End Results (SEER) Program
#' (\url{https://seer.cancer.gov/}) SEER*Stat Database.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/pharmacoevidence/seerSurv}
#'   \item \url{https://github.com/pharmacoevidence/seerSurv/issues}
#' }
#'
#' @author
#' Sameer Mansoori, Shubhram Pandey, Rashi Rani, Barinder Singh, Murat Kurt
#'
#' Maintainer: Shubhram Pandey \email{shubhram.pandey@@pharmacoevidence.com}
#'
#' @keywords internal
"_PACKAGE"
