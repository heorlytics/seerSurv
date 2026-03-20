#' US-CDC Life-Table for SEER-17 Background Mortality Adjustment
#'
#' @description
#' Annual death rates by single year of age (0–100) and sex from the US
#' Centers for Disease Control and Prevention (CDC) National Center for Health
#' Statistics, aligned to the SEER-17 study period (cases diagnosed before
#' 2021).  Used by \code{\link{make_background_surv}} and
#' \code{\link{run_tumour_analysis}} to normalise relative (disease-specific)
#' survival estimates to net (all-cause) survival.
#'
#' @format A data frame with 101 rows and 3 columns:
#' \describe{
#'   \item{\code{Age}}{Integer age (0, 1, 2, …, 100).}
#'   \item{\code{Males}}{Annual probability of death for males at this age.}
#'   \item{\code{Females}}{Annual probability of death for females at this age.}
#' }
#'
#' @details
#' A sex-mixed rate for a cohort with proportion male \eqn{\pi} is constructed
#' as:
#' \deqn{q_{\text{mix}} = \pi \cdot q_{\text{males}} + (1 - \pi) \cdot q_{\text{females}}.}
#' This is done inside \code{\link{run_tumour_analysis}} for each
#' recurrence-type subgroup.
#'
#' @source
#' US CDC/NCHS National Vital Statistics Reports, underlying cause-of-death
#' life-tables (\url{https://www.cdc.gov/nchs/nvss/life-tables.htm}).
#'
#' @seealso \code{\link{make_background_surv}}, \code{\link{run_tumour_analysis}}
#'
#' @examples
#' data(lifetable_seer)
#' head(lifetable_seer)
#' plot(lifetable_seer$Age, lifetable_seer$Males, type = "l",
#'      xlab = "Age", ylab = "Annual death rate", main = "US-CDC Life-Table")
#' lines(lifetable_seer$Age, lifetable_seer$Females, lty = 2)
#' legend("topleft", c("Males", "Females"), lty = 1:2)
"lifetable_seer"


#' SEER-17 Five-Year Aggregate Survival Data for 11 Cancer Types
#'
#' @description
#' Reference aggregate survival proportions at years 0–5 post-recurrence for
#' both local/regional-recurrence (LR) and distant-recurrence (DR) populations
#' across 11 cancer types, extracted from the Surveillance, Epidemiology, and
#' End Results SEER-17 registry (cases diagnosed before 2021).  Also includes
#' sample sizes, proportion male, and mean age for each subgroup.
#'
#' @format A data frame with 22 rows (2 recurrence types × 11 tumours) and
#' 12 columns:
#' \describe{
#'   \item{\code{Tumor}}{Cancer type label.}
#'   \item{\code{Recurrence}}{Recurrence type: \code{"LR"} or \code{"DR"}.}
#'   \item{\code{S0}}{Survival at year 0 (always 1).}
#'   \item{\code{S1}}{Survival at year 1.}
#'   \item{\code{S2}}{Survival at year 2.}
#'   \item{\code{S3}}{Survival at year 3.}
#'   \item{\code{S4}}{Survival at year 4.}
#'   \item{\code{S5}}{Survival at year 5.}
#'   \item{\code{prop_male}}{Proportion of patients who are male.}
#'   \item{\code{mean_age}}{Mean (or median) age at recurrence.}
#'   \item{\code{N_LR}}{Number of LR patients in the SEER cohort.}
#'   \item{\code{N_DR}}{Number of DR patients in the SEER cohort.}
#' }
#'
#' @details
#' The 11 tumour types are: Breast, Bladder, Colon & Rectum, Esophageal,
#' Kidney, Liver, Lung, Melanoma, Ovarian, Pancreatic, and Stomach.
#'
#' @source
#' Surveillance, Epidemiology, and End Results (SEER) Program
#' (\url{https://seer.cancer.gov/}) SEER*Stat Database:
#' SEER 17 Regs Research Data, Nov 2022 Sub (2000–2020).
#'
#' @seealso \code{\link{run_tumour_analysis}}
#'
#' @examples
#' data(tumour_data_seer)
#' dplyr::filter(tumour_data_seer, Tumor == "Melanoma")
"tumour_data_seer"
